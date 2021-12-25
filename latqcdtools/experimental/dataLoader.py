import numpy as np
import latqcdtools.logger as logger

class DataLoader:

    def __init__(self, path, use_sigma, fromNumpy=False,
             fromSample=False, sample_interval=None, standardize=False, precision='double'):
        self.path_train = path
        self.use_sigma = use_sigma
        self.fromNumpy = fromNumpy
        self.fromSample = fromSample
        self.sample_interval = sample_interval
        self.standardize = standardize
        
        if not self.use_sigma:
            self.sigma_y=None

        self.precision = np.float64 if precision == 'double' else np.float32

    
    def readFromNumpy(self):

        if self.fromSample:
            print(self.path_train)
            self.raw_samples = np.load(self.path_train)
            
            if self.raw_samples.shape[1] != 3:
                logger.warn("Data has unknown column format!")
                logger.warn("Full shape is:", self.raw_samples.shape)
            if self.use_sigma:
                train_x = self.raw_samples[:,-3]
                train_y = self.raw_samples[:,-2]
                sigma_y = self.raw_samples[:,-1]

            else:
                train_x = self.raw_samples[:,-3]
                train_y = self.raw_samples[:,-2]


            if self.sample_interval:
                train_x = train_x[self.sample_interval[0]:self.sample_interval[1]]
                train_y = train_y[self.sample_interval[0]:self.sample_interval[1]]
                if self.use_sigma:
                    sigma_y = sigma_y[self.sample_interval[0]:self.sample_interval[1]]
        else:
            print("Error: readFromNumpy with fromSample=False is not implemented yet.")
        
        if self.use_sigma:# or self.covar:
            return train_x, train_y, sigma_y
        return train_x, train_y


    def readFromTxt(self):

        if self.fromSample:
            print("Error: readFromTxt with fromSample=True is not implemented yet.")
        else:

            data = np.loadtxt(self.path_train).T
            if data.shape[0] != 3:
                logger.warn("Data has unknown column format!")
            train_x = data[0]
            train_y = data[1]

            if self.use_sigma:
                sigma_y = data[2]
            
        
        if self.use_sigma:
            return train_x, train_y, sigma_y
        return train_x, train_y

    

    def load(self):

        if self.fromNumpy:
            data = self.readFromNumpy()
        else:
            data = self.readFromTxt()
        
        train_x = data[0]
        train_y = data[1]
        if self.use_sigma:# or self.covar:
            sigma_y = data[2]

        if self.fromSample:
            dim = 1
        else:
            dim = 0

        
        means_x = train_x.mean(axis=dim, keepdims=True)
        stds_x = train_x.std(axis=dim, keepdims=True) + 1e-14 # prevent dividing by 0

        means_y = train_y.mean(axis=dim, keepdims=True)
        stds_y = train_y.std(axis=dim, keepdims=True) + 1e-14 # prevent dividing by 0
            
        if self.standardize:
            train_x = (train_x - means_x) / stds_x
            train_y = (train_y - means_y) / stds_y
            if self.use_sigma:
                sigma_y = sigma_y / stds_y
        else:
            means_x = np.zeros(means_x.shape)
            stds_x = np.ones(means_x.shape)

            means_y = np.zeros(means_y.shape)
            stds_y = np.ones(means_y.shape)
       

        self.train_x = train_x 
        self.train_y = train_y 
        if self.use_sigma:# or self.covar:
            self.sigma_y = sigma_y 

        self.means_x = means_x
        self.stds_x = stds_x 
       
        self.means_y = means_y
        self.stds_y = stds_y 

        

    def postProcess(self, x, y, sigma=[]):

        x = (x * self.stds_x + self.means_x)
        y = (y * self.stds_y + self.means_y)
        
        if len(sigma) != 0:
            sigma = sigma * self.stds_y
            return x, y, sigma
        return x, y
    
    def getTrainData(self):
        if self.use_sigma: # or self.covar:
            return self.postProcess(self.train_x,self.train_y,self.sigma_y)
        else:
            return self.postProcess(self.train_x,self.train_y)

    def createPredictData(self, test_x):
        
        if self.fromSample:
            test_x = np.repeat(test_x[np.newaxis],self.train_x.shape[0],axis=0)[:,:,np.newaxis]
        
        test_x = (test_x - self.means_x) / self.stds_x
        return test_x
