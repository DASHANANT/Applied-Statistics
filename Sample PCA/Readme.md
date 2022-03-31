# Principal Component Analysis based on Sample Covariance matrix

# Converting a Covariance Matrix to a Correlation Matrix

You can use similar operations to convert a covariance matrix to a correlation matrix. First, use the DIAG function to extract the variances from the diagonal elements of the covariance matrix. Then invert the matrix to form the diagonal matrix with diagonal elements that are the reciprocals of the standard deviations.
/** convert covariance matrix to correlation matrix **/
S = {1.0  1.0  8.1,
     1.0 16.0 18.0,
     8.1 18.0 81.0 };
 
/** standard deviations of each variable **/
D = sqrt(diag(S));
DInv = inv(D);
R = DInv * S * Dinv; /** correlation matrix **/
print R;
