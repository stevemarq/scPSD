import numpy as np
from scipy.fft import fftn

class scPSD:
    def __init__(self, scOmics):
        self.omics = np.asarray(scOmics)
        self.num_features, self.num_samples = self.omics.shape


    def correlation_est(self, data):
        covariance = np.cov(data)
        std = np.std(data, axis=1)[:, np.newaxis]
        std_mat = np.dot(std, std.T)
        p = covariance / std_mat
        a1 = p @ data
        return a1


    def genomic_fft(self, data):
        dft = fftn(np.abs(data), axes=0)
        magnitude = np.abs(dft)
        n = magnitude.shape[0]
        a2 = magnitude / n
        return a2


    def spectral_entropy_est(self, data):
        feature_sums = np.sum(data, axis=0)[:, np.newaxis]
        p_kj = (data.T * feature_sums[:, None]).T
        i_kj = -np.log(p_kj)
        h = p_kj * i_kj
        sample_means = np.mean(h, axis=1)[:, np.newaxis]
        a3 = -h + sample_means
        norm_a3 = (a3 - np.min(a3)) / (np.max(a3) - np.min(a3))

        return norm_a3


    def transform(self):
        correlation_est = self.correlation_est(self.omics)
        genomic_fft = self.genomic_fft(correlation_est)
        spectral_entropy_est = self.spectral_entropy_est(genomic_fft)
        return spectral_entropy_est

