import numpy as np
import scipy.integrate as integrate
import scipy.optimize
from matplotlib import pyplot as plt

def gauss(x, params):
    """
    :param x: input data
    :param params: list of parameters for the model
    :return:
    """
    p = params
    result = p[0] * np.exp(-np.power(x - p[1], 2) / p[2])
    return result


def linear(x, params):
    """
    :param x: input data
    :param params: list of parameters for the model
    :return:
    """
    p = params
    result = p[0] + x * p[1]
    return result


def lorrentz(x, params):
    """
    :param x: input data
    :param params: list of parameters for the model
    :return:
    """

    g = params[2]
    x0 = params[1]
    scale = params[0]

    return scale * g / (np.power(x - x0, 2) + np.power(g, 2))

class fit:
    """
    Class to handle curve fitting with arbitrary curve shapes and peak integration
    User initialize class with a pair of [x,y] where x and y are the
    xl and xu set the lower and upper limits for curve fitting based on index in the input data
    Usage:
        use addmodel() to add desired models by sucessively calling
        use addbckg() to add background
        use buildmodel() to construct the model for fitting
        use fit() to fit the model with input be a list of parameters for each model in the order it was added
        use plot() to plot the resulting models
        use getparams() to get the resulting parameters from the fit
        use integrate() to integrate the curves from the fit

    If addmodel() or addbckg() is called, the model must be rebuilt and fit again
    """
    def __init__(self, pair, xl=200, xu=600):
        self.x = pair[0]
        self.y = pair[1]
        self.params = None
        self.size = len(self.x)
        self.restrict = None
        self.model = []
        self.modelinitial = None
        self.bckg = None
        self.bckgparams = None
        self.full = None
        self.modelparams = None
        self.isBuilt = False
        self.isFit = False
        self.xl = xl
        self.xu = xu
        self.modelparamcount = 0
        self.integrated = []


    def addmodel(self, fn, params):
        """
        Add a model that will be linearly added to the other models and background
        :param fn: a function that has as input a,b where a is a 1D vector of data to predict
        and b is a 1D vector of parameters for the model
        :param params: the number of parameters that the function takes
        :return: None
        """
        if len(self.model) == 0:
            self.model = [fn]
            self.modelparams = [params]
            self.isBuilt = False
        else:
            self.model.append(fn)
            self.modelparams.append(params)
            self.isBuilt = False

    def addbckg(self, shape, params):
        """
        Add a background that will be linearly added to the models
        :param fn: a function that has as input a,b where a is a 1D vector of data to predict
        and b is a 1D vector of parameters for the model
        :param params: the number of parameters that the function takes
        :return: None
        """
        self.bckg = shape
        self.bckgparams = params
        self.isBuilt = False


    def buildmodel(self):
        """
        Assembles the full model by adding the models and backgrounds
        :return: None
        """
        if self.model is None:
            raise ValueError("model is not set")
        if self.bckg is None:
            raise ValueError("background is not set")

        def f(x, *params):
            parampos = 0
            result = self.model[0](x, params[parampos:self.modelparams[0]])
            parampos += self.modelparams[0]

            for i in range(1, len(self.model)):
                result += self.model[i](x, params[parampos:parampos + self.modelparams[i]])
                parampos += self.modelparams[i]

            result += self.bckg(x, params[parampos:])

            return result  # self.model(x, params[:self.modelparams]) + self.bckg(x, params[self.modelparams:])

        self.modelparamcount = sum(self.modelparams)

        self.full = f
        self.isBuilt = True
        self.isFit = False

    def getparams(self):
        """
        :return: The resulting parameters from fitting as a 2D numpy array of form [data index, parameter index]
        """
        if not (self.isFit and self.isBuilt):
            raise AttributeError("Current model has not been fit")
        return self.params

    def fit(self, guess, bguess):
        """
        Fit the full model
        :param guess: parameters for the models to initialize using
        Parameters are a list in corresponding to the same order the models were added in
        :param bguess: parameters to initialize the background to
        :return: None
        """
        if not self.isBuilt:
            raise AttributeError("model has not yet been built")
        if len(guess) != self.modelparamcount:
            raise ValueError("wrong model param dim")
        if len(bguess) != self.bckgparams:
            raise ValueError("wrong bckg param dim")

        self.params = []

        for i in range(self.size):

            fit = scipy.optimize.curve_fit(self.full, self.x[i][self.xl:self.xu], self.y[i][self.xl:self.xu],
                                           p0=guess + bguess)
            self.params.append(fit[0])
        print("fit " + str(self.size) + " models")
        self.isFit = True

    def getmodel(self, i, modelnumber, bckgsign=0):
        """
        Returns a fit model (and background if desired) for evaluation on data
        :param i: the index of data to get the model for
        :param modelnumber: the zero indexed model to return
        :param bckgsign: multiples the background by sign. If 0, does not include the background
        :return: a function that evaluates a model (can include background based on bckg sign) using the
        parameters from the fitting
        """
        def individualmodel(x):
            bckg = self.bckg(x, self.getparams()[i][self.modelparamcount:])
            paramstart = sum(self.modelparams[:modelnumber])
            paramend = sum(self.modelparams[:modelnumber + 1])
            # print(paramstart, paramend)
            model = self.model[modelnumber](x, self.getparams()[i][paramstart:paramend])
            return model + (bckg * bckgsign if bckgsign else 0)

        return individualmodel

    def plot(self):
        """
        For each dataset make a subplot with the original data, the fit model, the backgtround, and inidividual curves
        :return: None
        """
        if not (self.isFit and self.isBuilt):
            raise AttributeError("Current model has not been fit")

        fig = plt.figure(figsize=(10, 60))
        for i in range(self.size):
            ax = fig.add_subplot(23, 3, i + 1)
            ax.plot(self.x[i][self.xl:self.xu], self.y[i][self.xl:self.xu])
            ax.plot(self.x[i][self.xl:self.xu], self.full(self.x[i][self.xl:self.xu], *self.getparams()[i]))
            ax.plot(self.x[i][self.xl:self.xu],
                    self.bckg(self.x[i][self.xl:self.xu], self.getparams()[i][self.modelparamcount:]))


            def plotcurves(self, x, idx, ax):

                for modelidx in range(len(self.model)):
                    model = self.getmodel(idx, modelidx, 1)
                    ax.plot(x[i][self.xl:self.xu], model(x[i][self.xl:self.xu]))

            plotcurves(self, self.x, i, ax)

            ax.text(0.5, .9, i, transform=ax.transAxes)
        plt.show()



    def integrate(self):
        """
        Integrates the fit models using trapezoild integration
        :return: Returns a numpy array of form [data index, model number]
        """
        self.integrated = np.zeros([len(self.model), len(self.x)])
        for i in range(len(self.x)):
            parampos = 0
            x = self.x[i]
            for m in range(len(self.model)):
                integrand = self.getmodel(i, m, bckgsign=0)
                y = integrand(x)
                self.integrated[m, i] = integrate.cumtrapz(y, x)[-1]

        return self.integrated

