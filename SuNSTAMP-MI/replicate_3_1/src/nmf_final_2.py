import numbers
import os
import pandas as pd
import time
from numpy import sqrt
import numpy as np
import scipy.sparse as sp
from scipy import linalg, sparse
from scipy import linalg



# Set the seed for NumPy so that each job id produces same input/mask sequences
np.random.seed(42)
task_id = int((os.getenv('SLURM_ARRAY_TASK_ID')))


# Global mask matrix
global_mask_matrix = None
global_reverse_mask_matrix = None
# Global Parameters
EPSILON = np.finfo(np.float32).eps
LAMBDA= 0.0000275

months = 3
replicate = 1


def safe_sparse_dot(a, b, *, dense_output=False):
    """Dot product that handle the sparse matrix case correctly.

    Parameters
    ----------
    a : {ndarray, sparse matrix}
    b : {ndarray, sparse matrix}
    dense_output : bool, default=False
        When False, ``a`` and ``b`` both being sparse will yield sparse output.
        When True, output will always be a dense array.

    Returns
    -------
    dot_product : {ndarray, sparse matrix}
        Sparse if ``a`` and ``b`` are sparse and ``dense_output=False``.

    Examples
    --------
    >>> from scipy.sparse import csr_matrix
    >>> from sklearn.utils.extmath import safe_sparse_dot
    >>> X = csr_matrix([[1, 2], [3, 4], [5, 6]])
    >>> dot_product = safe_sparse_dot(X, X.T)
    >>> dot_product.toarray()
    array([[ 5, 11, 17],
           [11, 25, 39],
           [17, 39, 61]])
    """
    if a.ndim > 2 or b.ndim > 2:
        if sparse.issparse(a):
            # sparse is always 2D. Implies b is 3D+
            # [i, j] @ [k, ..., l, m, n] -> [i, k, ..., l, n]
            b_ = np.rollaxis(b, -2)
            b_2d = b_.reshape((b.shape[-2], -1))
            ret = a @ b_2d
            ret = ret.reshape(a.shape[0], *b_.shape[1:])
        elif sparse.issparse(b):
            # sparse is always 2D. Implies a is 3D+
            # [k, ..., l, m] @ [i, j] -> [k, ..., l, j]
            a_2d = a.reshape(-1, a.shape[-1])
            ret = a_2d @ b
            ret = ret.reshape(*a.shape[:-1], b.shape[1])
        else:
            ret = np.dot(a, b)
    else:
        ret = a @ b

    if (
        sparse.issparse(a)
        and sparse.issparse(b)
        and dense_output
        and hasattr(ret, "toarray")
    ):
        return ret.toarray()
    return ret

def squared_norm(x):
    """Squared Euclidean or Frobenius norm of x.

    Faster than norm(x) ** 2.

    Parameters
    ----------
    x : array-like
        The input array which could be either be a vector or a 2 dimensional array.

    Returns
    -------
    float
        The Euclidean norm when x is a vector, the Frobenius norm when x
        is a matrix (2-d array).
    """
    x = np.ravel(x, order="K")
    if np.issubdtype(x.dtype, np.integer):
        print("Array type is integer, np.dot may overflow. Data should be float type to avoid this issue")
          
    return np.dot(x, x)

def norm(x):
    """Dot product-based Euclidean norm implementation.

    See: http://fa.bianp.net/blog/2011/computing-the-vector-norm/

    Parameters
    ----------
    x : array-like
        Vector for which to compute the norm.
    """
    return sqrt(squared_norm(x))


def trace_dot(X, Y):
    """Trace of np.dot(X, Y.T).

    Parameters
    ----------
    X : array-like
        First matrix.
    Y : array-like
        Second matrix.
    """
    return np.dot(X.ravel(), Y.ravel())

def _beta_loss_to_float(beta_loss):
    """Convert string beta_loss to float."""
    beta_loss_map = {"frobenius": 2, "kullback-leibler": 1, "itakura-saito": 0}
    if isinstance(beta_loss, str):
        beta_loss = beta_loss_map[beta_loss]
    return beta_loss

def _beta_divergence(X, W, H, beta, square_root=False):
    # Frobenius norm
    updated_loss = global_mask_matrix * (X - np.dot(W,H))
    res = squared_norm(updated_loss)
    return res

def frobenius_norm_first_column(W, Ideal_Impact_Matrix):
    """
    Compute the Frobenius norm of the difference between the first column of W and the Ideal_Impact_Matrix.

    Parameters:
    - W (np.ndarray): A matrix of dimensions (n, r).
    - Ideal_Impact_Matrix (np.ndarray): A matrix of dimensions (n, 1).

    Returns:
    - float: The Frobenius norm of the difference.
    """
    # Compute the difference between the first column of W and Ideal_Impact_Matrix
    difference = LAMBDA*(Ideal_Impact_Matrix.ravel()) - W[:, 0]

    # Calculate the Frobenius norm of the difference
    updated_difference = (global_mask_matrix[:,0]) * difference
    norm = squared_norm(updated_difference) 
    return norm
    

def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance.

    Parameters
    ----------
    seed : None, int or instance of RandomState
        If seed is None, return the RandomState singleton used by np.random.
        If seed is an int, return a new RandomState instance seeded with seed.
        If seed is already a RandomState instance, return it.
        Otherwise raise ValueError.

    Returns
    -------
    :class:`numpy:numpy.random.RandomState`
        The random state object based on `seed` parameter.

    Examples
    --------
    >>> from sklearn.utils.validation import check_random_state
    >>> check_random_state(42)
    RandomState(MT19937) at 0x...
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, numbers.Integral):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError(
        "%r cannot be used to seed a numpy.random.RandomState instance" % seed
    )


    
def _initialize_nmf(X, n_components, init=None, eps=1e-6, random_state=None):
    """Algorithms for NMF initialization.

    Computes an initial guess for the non-negative
    rank k matrix approximation for X: X = WH.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        The data matrix to be decomposed.

    n_components : int
        The number of components desired in the approximation.

    init :  {'random', 'nndsvd', 'nndsvda', 'nndsvdar'}, default=None
        Method used to initialize the procedure.
        Valid options:

        - None: 'nndsvda' if n_components <= min(n_samples, n_features),
            otherwise 'random'.

        - 'random': non-negative random matrices, scaled with:
            sqrt(X.mean() / n_components)


    eps : float, default=1e-6
        Truncate all values less then this in output to zero.

    random_state : int, RandomState instance or None, default=None
        Used when ``init`` == 'nndsvdar' or 'random'. Pass an int for
        reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    W : array-like of shape (n_samples, n_components)
        Initial guesses for solving X ~= WH.

    H : array-like of shape (n_components, n_features)
        Initial guesses for solving X ~= WH.

    References
    ----------
    C. Boutsidis, E. Gallopoulos: SVD based initialization: A head start for
    nonnegative matrix factorization - Pattern Recognition, 2008
    http://tinyurl.com/nndsvd
    """
    n_samples, n_features = X.shape
    # Random initialization
    if init == "random":
        avg = np.sqrt(X.mean() / n_components)
        rng = check_random_state(random_state)
        H = avg * rng.standard_normal(size=(n_components, n_features)).astype(
            X.dtype, copy=False
        )
        W = avg * rng.standard_normal(size=(n_samples, n_components)).astype(
            X.dtype, copy=False
        )
        np.abs(H, out=H)
        np.abs(W, out=W)
        return W, H
    


def _multiplicative_update_w(
    X,
    W,
    H,
    target_values,
    beta_loss,
    l1_reg_W,
    l2_reg_W,
    gamma,
    H_sum=None,
    HHt=None,
    XHt=None,
    update_H=True,
):
    """Update W in Multiplicative Update NMF."""
    if beta_loss == 2:
        # Numerator
        if XHt is None:
            XHt = safe_sparse_dot(X, H.T)
        if update_H:
            # avoid a copy of XHt, which will be re-computed (update_H=True)
            numerator = XHt
        else:
            # preserve the XHt, which is not re-computed (update_H=False)
            numerator = XHt.copy()

        # Denominator
        if HHt is None:
            HHt = np.dot(H, H.T)
        denominator = np.dot(W, HHt)

	# ELSE WILL NEVER APPEAR

    # Add L1 and L2 regularization
    if l1_reg_W > 0:
        denominator += l1_reg_W
    if l2_reg_W > 0:
        denominator = denominator + l2_reg_W * W
    denominator[denominator == 0] = EPSILON

    

    # -------------------------- MAKING THE FIRST COLUMN TO BE ZERO ---------
    # Set the first column to all zeros
    #delta_W[:, 0] = 0.5 * (delta_W[:, 0])
    # -------------------------- MAKING THE FIRST COLUMN TO BE ZERO ---------
    
	# --------------------- ADDING SUPERVISION HERE -------------
 
    # Supervision: Initialize delta_W_2 to zero matrix
    delta_W_2 = np.zeros_like(W)
        
    # Supervise the first column: division of target_values by the first column of W
    supervised_update = target_values
    #supervised_update[np.isnan(supervised_update)] = 0  # Handling possible NaNs
        
    # Assign the supervised update to the first column of delta_W_2
    delta_W_2[:, 0] = supervised_update
    delta_W_2 = delta_W_2 * LAMBDA

    # Add the supervised update delta_W_2 to delta_W
    #delta_W += delta_W_2
	# --------------------- ADDING SUPERVISION HERE -------------


    numerator += (delta_W_2)
    denominator += (W)
    numerator /= denominator
    delta_W = numerator

    # gamma is in ]0, 1]
    if gamma != 1:
        delta_W **= gamma

    W *= delta_W

    return W, H_sum, HHt, XHt


def _multiplicative_update_h(
    X, W, H, beta_loss, l1_reg_H, l2_reg_H, gamma, A=None, B=None, rho=None
):
    """update H in Multiplicative Update NMF."""
    if beta_loss == 2:
        numerator = safe_sparse_dot(W.T, X)
        denominator = np.linalg.multi_dot([W.T, W, H])

    # ELSE WILL NEVER APPEAR

    # Add L1 and L2 regularization
    if l1_reg_H > 0:
        denominator += l1_reg_H
    if l2_reg_H > 0:
        denominator = denominator + l2_reg_H * H
    denominator[denominator == 0] = EPSILON

    if A is not None and B is not None:
        # Updates for the online nmf
        if gamma != 1:
            H **= 1 / gamma
        numerator *= H
        A *= rho
        B *= rho
        A += numerator
        B += denominator
        H = A / B

        if gamma != 1:
            H **= gamma
    else:
        delta_H = numerator
        delta_H /= denominator
        if gamma != 1:
            delta_H **= gamma
        H *= delta_H

    return H


def _fit_multiplicative_update(
    X,
    W,
    H,
    target_values,
    beta_loss="frobenius",
    max_iter=200,
    tol=1e-3,
    l1_reg_W=0,
    l1_reg_H=0,
    l2_reg_W=0,
    l2_reg_H=0,
    update_H=True,
    verbose=1,
):
    """Compute Non-negative Matrix Factorization with Multiplicative Update.

    The objective function is _beta_divergence(X, WH) and is minimized with an
    alternating minimization of W and H. Each minimization is done with a
    Multiplicative Update.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Constant input matrix.

    W : array-like of shape (n_samples, n_components)
        Initial guess for the solution.

    H : array-like of shape (n_components, n_features)
        Initial guess for the solution.
    target_values: the ideal impact values/ attraction values

    beta_loss : float or {'frobenius', 'kullback-leibler', \
            'itakura-saito'}, default='frobenius'
        String must be in {'frobenius', 'kullback-leibler', 'itakura-saito'}.
        Beta divergence to be minimized, measuring the distance between X
        and the dot product WH. Note that values different from 'frobenius'
        (or 2) and 'kullback-leibler' (or 1) lead to significantly slower
        fits. Note that for beta_loss <= 0 (or 'itakura-saito'), the input
        matrix X cannot contain zeros.

    max_iter : int, default=200
        Number of iterations.

    tol : float, default=1e-4
        Tolerance of the stopping condition.

    l1_reg_W : float, default=0.
        L1 regularization parameter for W.

    l1_reg_H : float, default=0.
        L1 regularization parameter for H.

    l2_reg_W : float, default=0.
        L2 regularization parameter for W.

    l2_reg_H : float, default=0.
        L2 regularization parameter for H.

    update_H : bool, default=True
        Set to True, both W and H will be estimated from initial guesses.
        Set to False, only W will be estimated.

    verbose : int, default=0
        The verbosity level.

    Returns
    -------
    W : ndarray of shape (n_samples, n_components)
        Solution to the non-negative least squares problem.

    H : ndarray of shape (n_components, n_features)
        Solution to the non-negative least squares problem.

    n_iter : int
        The number of iterations done by the algorithm.

    References
    ----------
    Lee, D. D., & Seung, H., S. (2001). Algorithms for Non-negative Matrix
    Factorization. Adv. Neural Inform. Process. Syst.. 13.
    Fevotte, C., & Idier, J. (2011). Algorithms for nonnegative matrix
    factorization with the beta-divergence. Neural Computation, 23(9).
    """
    start_time = time.time()

    beta_loss = _beta_loss_to_float(beta_loss)

    # gamma for Maximization-Minimization (MM) algorithm [Fevotte 2011]
    if beta_loss < 1:
        gamma = 1.0 / (2.0 - beta_loss)
    elif beta_loss > 2:
        gamma = 1.0 / (beta_loss - 1.0)
    else:
        gamma = 1.0

    # used for the convergence criterion
    error_at_init = _beta_divergence(X, W, H, beta_loss, square_root=True)
    error_at_init_first_column = frobenius_norm_first_column(W, target_values)
    previous_error = error_at_init + error_at_init_first_column

    H_sum, HHt, XHt = None, None, None
    for n_iter in range(1, max_iter + 1):
        # update W
        # H_sum, HHt and XHt are saved and reused if not update_H
        W, H_sum, HHt, XHt = _multiplicative_update_w(
            X,
            W,
            H,
            target_values=target_values,
            beta_loss=beta_loss,
            l1_reg_W=l1_reg_W,
            l2_reg_W=l2_reg_W,
            gamma=gamma,
            H_sum=H_sum,
            HHt=HHt,
            XHt=XHt,
            update_H=update_H,
        )

        # necessary for stability with beta_loss < 1
        if beta_loss < 1:
            W[W < np.finfo(np.float64).eps] = 0.0

        # update H (only at fit or fit_transform)
        if update_H:
            H = _multiplicative_update_h(
                X,
                W,
                H,
                beta_loss=beta_loss,
                l1_reg_H=l1_reg_H,
                l2_reg_H=l2_reg_H,
                gamma=gamma,
            )

            # These values will be recomputed since H changed
            H_sum, HHt, XHt = None, None, None

            # necessary for stability with beta_loss < 1
            if beta_loss <= 1:
                H[H < np.finfo(np.float64).eps] = 0.0

        # test convergence criterion every 1 iterations
        if tol > 0 and n_iter % 1 == 0:
            error_beta_divergence = _beta_divergence(X, W, H, beta_loss, square_root=True)
            error_first_col = frobenius_norm_first_column(W, target_values)
            error = error_beta_divergence + error_first_col

            if verbose:
                iter_time = time.time()
                #print(
                    #"Epoch %02d reached after %.3f seconds, error: %f"
                    #% (n_iter, iter_time - start_time, error)
                #)
            #print(previous_error)
            #print(error)
            if (previous_error - error) / error_at_init < tol:
                break
            previous_error = error

    # do not print if we have already printed in the convergence test
    if verbose and (tol == 0 or n_iter % 10 != 0):
        end_time = time.time()
        #print(
            #"Epoch %02d reached after %.3f seconds." % (n_iter, end_time - start_time)
        #)
    
    #print(W)
    #print("W matrix printed")
    #print(H)
    #print("H matrix printed")
    
    return W, H, n_iter



# Load gene expression data
# gene_df = pd.read_csv('gene_exp.csv', index_col=0)
gene_df = pd.read_csv(f'../data/gene_exp_{months}_{replicate}.csv', index_col=0)

# Append the data to the CSV file



# Load target values (total attraction values) which also contains coordinates
# coord_df = pd.read_csv('normalized_cell_attraction_results.csv')
coord_df = pd.read_csv(f'../data/supervision_{months}_{replicate}.csv')

# Perform an inner join to keep only cells present in both DataFrames
common_df = gene_df.join(coord_df.set_index('Cell_Id'), how='inner')

# If the coordinates are required separately you can extract them, assuming 'X' and 'Y' are columns
coordinates = common_df[['X-scaled', 'Y-scaled']]


# Load and preprocess target values (total attraction values)
target_values = common_df['Total_Attraction']

# Remove coordinate columns if not needed in gene expression data
gene_df = common_df.drop(['X-scaled', 'Y-scaled','Total_Attraction','Normalized_Total_Attraction'], axis=1)



# Iterate over the number of components from 5 to 20
data_matrix = gene_df.to_numpy()

gene_names = gene_df.columns[0:] 



############################## Nested k-fold cross validation ###################

# Define the number folds or test sets
K = 10
# lambdas = [3e-6, 2.75e-6, 2.5e-6, 2e-6, 1.75e-6, 1.5e-6, 1.25e-6, 1e-6,
#            3e-5, 2.75e-5, 2.5e-5, 2e-5, 1.75e-5, 1.5e-5, 1.25e-5, 1e-5,
#            3e-4, 2.75e-4, 2.5e-4, 2e-4, 1.75e-4, 1.5e-4, 1.25e-4, 1e-4,
#            1e-3, 1e-2]
lambdas = [
    0.01, 0.02, 0.05,
    0.1, 0.2, 0.5,
    1, 2, 5,
    8, 10, 12,
    14, 15, 16,
    18, 20, 25,
    30, 40, 50,
    70, 100, 200,
    500, 1000
]
# Define ranges or sets of values for each hyperparameter
# rank_values = []
# if (task_id == 0):
#     rank_values=[11]
# elif (task_id==1):
#     rank_values=[12]
# elif (task_id==2):
#     rank_values=[13]
# elif (task_id==3):
#     rank_values=[14]
# elif (task_id==4):
#     rank_values=[15]
# elif (task_id==5):
#     rank_values=[16]
# elif (task_id==6):
#     rank_values=[17]
# elif (task_id==7):
#     rank_values=[18]
# elif (task_id==8):
#     rank_values=[19]
# elif (task_id==9):
#     rank_values=[20]
# elif (task_id==10):
#     rank_values=[21]
# elif (task_id==11):
#     rank_values=[22]
# elif (task_id==12):
#     rank_values=[23]
# elif (task_id==13):
#     rank_values=[24]
# elif (task_id==14):
#     rank_values=[25]
rank_values = []
if task_id == 0:
    rank_values = [5]
elif task_id == 1:
    rank_values = [6]
elif task_id == 2:
    rank_values = [7]
elif task_id == 3:
    rank_values = [8]
elif task_id == 4:
    rank_values = [9]
elif task_id == 5:
    rank_values = [10]
elif task_id == 6:
    rank_values = [11]
elif task_id == 7:
    rank_values = [12]
elif task_id == 8:
    rank_values = [13]
elif task_id == 9:
    rank_values = [14]
elif task_id == 10:
    rank_values = [15]
elif task_id == 11:
    rank_values = [16]
elif task_id == 12:
    rank_values = [17]
elif task_id == 13:
    rank_values = [18]
elif task_id == 14:
    rank_values = [19]
elif task_id == 15:
    rank_values = [20]
elif task_id == 16:
    rank_values = [21]
elif task_id == 17:
    rank_values = [22]
elif task_id == 18:
    rank_values = [23]
elif task_id == 19:
    rank_values = [24]
elif task_id == 20:
    rank_values = [25]

    

data_list = []
final_loss=float('inf')
final_rank = 5
final_lambda = 1
# Create a dictionary to store the mapping from hyperparameters to loss values
#hyperparameter_losses = {(p1, p2): [] for p1, p2 in product(rank_values, lambdas)}

hyperparameter_losses = {}
for p1 in rank_values:
    for p2 in lambdas:
        hyperparameter_losses[(p1, p2)] = []
        

mask_matrices_outer = np.ones((K, *data_matrix.shape))


num_rows = data_matrix.shape[0]
test_set_size = int(0.1 * num_rows)

all_omitted_indices = set()

for i in range(K):
    
    mask_matrix_outer = np.ones(data_matrix.shape)
    
    
    available_indices = list(set(range(num_rows)) - all_omitted_indices)
    
    # If there are not enough available indices, reset the omitted set
    if len(available_indices) < test_set_size:
        all_omitted_indices.clear()
        available_indices = list(set(range(num_rows)))
    
    # Randomly choosing indices for the test set
    mask_indices_outer = np.random.choice(available_indices, size=test_set_size, replace=False)
    
    # Addding the chosen indices to the omitted set
    all_omitted_indices.update(mask_indices_outer)
    
    # chosen indices to zero to omit them
    mask_matrix_outer[mask_indices_outer, :] = 0
    
    # Storing the outer mask matrix in the array
    mask_matrices_outer[i] = mask_matrix_outer


final_loss_list = []
start = 1
for outer_mask in mask_matrices_outer:
	# -------------- reverse outer mask -----------
    reverse_outer_mask = np.zeros(data_matrix.shape)
    reverse_outer_indices = np.where(outer_mask[:, 0] == 0)[0]
    reverse_outer_mask[reverse_outer_indices, :] = 1

	# -------------- reverse outer mask -----------
    # rows that are not omitted in the outer mask
    remaining_indices = np.where(outer_mask[:, 0] == 1)[0]
    num_remaining_rows = len(remaining_indices)
    inner_test_set_size = int(0.1 * num_remaining_rows)

    # Initialize an array to store the inner mask matrices for the current outer mask
    inner_masks_for_current_outer = np.ones((K, *data_matrix.shape))

	
    # Keep track of already omitted indices for the inner loop
    all_inner_omitted_indices = set()
    mask_matrices_inner = []
    all_losses=[]
	
    # Generate the inner mask matrices
    for j in range(K): 
        # Start with the outer mask matrix
        # rows that are already omitted in outer mask will also be omitted here
        mask_matrix_inner = np.copy(outer_mask)
        
        # Find available indices that have not been omitted yet in the inner loop
        available_inner_indices = list(set(remaining_indices) - all_inner_omitted_indices)
        
        # If there are not enough available indices, reset the omitted set
        if len(available_inner_indices) < inner_test_set_size:
            all_inner_omitted_indices.clear()
            available_inner_indices = list(remaining_indices)
        
        # Randomly choose indices for the inner test set
        mask_indices_inner = np.random.choice(available_inner_indices, size=inner_test_set_size, replace=False)
        
        # Add the chosen indices to the omitted set
        all_inner_omitted_indices.update(mask_indices_inner)
        
        # Set the chosen indices to zero to omit them
        mask_matrix_inner[mask_indices_inner, :] = 0
        
        # Store the inner mask matrix in the array for the current outer mask
        inner_masks_for_current_outer[j] = mask_matrix_inner
        # Append the inner masks for the current outer mask to the main list
        mask_matrices_inner.append(inner_masks_for_current_outer[j])
        
	
	# the loop will be run here
    for inner_mask in mask_matrices_inner:
        # do jobs here
        reverse_inner_mask = np.zeros(data_matrix.shape)
        reverse_indices = np.where(inner_mask[:, 0] == 0)[0]
        reverse_inner_mask[reverse_indices, :] = 1
        # Reverse Inner Matrix Formed
        for num_components in rank_values:
            for lambda_now in lambdas:
				# Initialize and fit the NMF model
                global_mask_matrix = inner_mask
                #print(global_mask_matrix)
                
                LAMBDA = lambda_now
                W, H = _initialize_nmf(data_matrix, n_components=num_components, init='random')
                W, H, n_iter = _fit_multiplicative_update(data_matrix, W=W, H=H, target_values=target_values)
                

				# Create a DataFrame from the W matrix
                w_df = pd.DataFrame(W, columns=[f'Component_{i+1}' for i in range(num_components)])
                w_df.insert(0, 'Cell_Id', common_df.index)
                # Save the DataFrame to a CSV file
                #w_df.to_csv(f'output_W_matrix_rank_{num_components}.csv', index=False)
                #print(f'Output W matrix for rank {num_components} saved as output_W_matrix_rank_{num_components}.csv')
                

                # ----------------------- OUTPUT THE 1ST ROW OF H MATRIX-----------------
                h_first_row = H[0, :]
                h_df = pd.DataFrame({'Gene': gene_names, 'Value': h_first_row})
                # Sort the DataFrame by the 'Value' column
                # Save the H DataFrame to a CSV file
                #h_df.to_csv(f'output_H_matrix_first_row_rank_{num_components}.csv', index=False)
                
				# ---------- TESTING LOSS ------------------------
                now_loss = reverse_inner_mask * (data_matrix - np.dot(W,H))
                frob_norm_data = np.linalg.norm(now_loss, ord='fro')
                # Append the loss to the list in the dictionary
                all_losses.append((num_components, lambda_now, frob_norm_data))
                hyperparameter_losses[(num_components, lambda_now)].append(frob_norm_data)
                

	# Outside the innner loop
    # Calculate the total loss for each combination of hyperparameters
    results_df = pd.DataFrame(all_losses, columns=['Rank', 'Lambda', 'Loss'])
    output_file_name_1 = f'../result/mou_errors_{task_id}_{months}_{replicate}.csv'
    results_df.to_csv(output_file_name_1, mode='a', header=False, index=False)
    all_losses=[]
    min_loss = float('inf')
    total_losses = {}
    best_hyperparameters = None
    for key, losses in hyperparameter_losses.items():
        total_losses[key] = sum(losses)
        #print(sum(losses))
        #print(losses)
        if min_loss > sum(losses):
            min_loss = sum(losses)
            best_hyperparameters = key
            
        
    
    optimal_rank_outer = best_hyperparameters[0]
    optimal_lambda_outer= best_hyperparameters[1]
    

	# ----------------	 MODEL TRAINIG OUTER LOOP ---------------
    # Initialize and fit the NMF model
    LAMBDA = optimal_lambda_outer
    global_mask_matrix = outer_mask
    W, H = _initialize_nmf(data_matrix, n_components=optimal_rank_outer, init='random')
    W, H, n_iter = _fit_multiplicative_update(data_matrix, W=W, H=H, target_values=target_values)
    # ---------- TESTING LOSS ------------------------
    now_loss = reverse_outer_mask * (data_matrix - np.dot(W,H))
    frob_norm_loss = np.linalg.norm(now_loss, ord='fro')
    final_loss_list.append((optimal_rank_outer, optimal_lambda_outer, frob_norm_loss))
    print(optimal_rank_outer, optimal_lambda_outer, frob_norm_loss)
    data_list.append((optimal_rank_outer, optimal_lambda_outer, frob_norm_loss))
    if final_loss > frob_norm_loss:
        final_loss = frob_norm_loss
        final_lambda = optimal_lambda_outer
        final_rank = optimal_rank_outer


# Create a pandas DataFrame from the list of tuples
df = pd.DataFrame(data_list, columns=["Optimal Rank Outer", "Optimal Lambda Outer", "Frobenius Norm Loss"])

# Write the DataFrame to a CSV file
output_file_name_2 = f'../result/mou_{task_id}_{months}_{replicate}.csv'
df.to_csv(output_file_name_2, index=False)

print("Data has been written to output.csv")

now_matrix = np.ones(data_matrix.shape)
global_mask_matrix = now_matrix
LAMBDA = final_lambda
W, H = _initialize_nmf(data_matrix, n_components=final_rank, init='random')
W, H, n_iter = _fit_multiplicative_update(data_matrix, W=W, H=H, target_values=target_values)
# Create a DataFrame from the W matrix
w_df = pd.DataFrame(W, columns=[f'Component_{i+1}' for i in range(final_rank)])
w_df.insert(0, 'Cell_Id', common_df.index)
# Save the DataFrame to a CSV file
w_df.to_csv(f'../result/output_W_matrix_rank_{final_rank}_{months}_{replicate}.csv', index=False)
print(f'Output W matrix for rank {final_rank} saved as output_W_matrix_rank_{final_rank}.csv')
    # ----------------------- OUTPUT THE 1ST ROW OF H MATRIX-----------------
    # Handling the H matrix for output
    # Extract the first row of the H matrix and create a DataFrame
h_first_row = H[0, :]
h_df = pd.DataFrame({'Gene': gene_names, 'Value': h_first_row})

    # Sort the DataFrame by the 'Value' column
h_df = h_df.sort_values(by='Value', ascending=False)

    # Save the H DataFrame to a CSV file
h_df.to_csv(f'../result/output_H_matrix_first_row_rank_{final_rank}_{months}_{replicate}.csv', index=False)
print(f'Output H matrix first row for rank {final_rank} sorted and saved as output_H_matrix_first_row_rank_{num_components}.csv')
    # TESTING LOSS FOR THOSE 205 DATA------------------------
now_loss = (data_matrix - np.dot(W,H))
frob_norm_data = np.linalg.norm(now_loss, ord='fro')
print("done: ")
print(frob_norm_data)



results_df = pd.DataFrame(data_list, columns=['Rank', 'Lambda','Loss'])
output_file_name_3 = f'../result/nmf_errors_mou_{task_id}_{months}_{replicate}.csv'
results_df.to_csv(output_file_name_3, index=False)

# https://machinelearningmastery.com/nested-cross-validation-for-machine-learning-with-python/