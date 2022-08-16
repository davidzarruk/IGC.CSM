
#' Array operator to mimic different-dimension-array element-wise operations in
#' MATLAB. It receives as input two arrays of potentially different dimensions,
#' it resizes them to have same dimensions and finally performs the element-wise
#' operation.
#'
#' @param array1 The first array
#' @param array2 The second array
#' @param operation The operation. It can take values: '+', '-', '*',
#' '/' and '^'
#'
#' @return An array with dimensions equal to the "largest" input array. It is
#' the result of applying the operator element-wise to both input arrays.
#'
#' @examples
array_operator = function(array1, array2, operation){
  if(is.array(array1) == FALSE){
    array1 = array(array1, dim=dim(array2))
  } else if(is.array(array2) == FALSE){
    array2 = array(array2, dim=dim(array1))
  }
  dim1 = dim(array1)
  dim2 = dim(array2)
  if(length(dim1) < length(dim2)){
    array1 = kronecker(array1, array(1, dim = array(1, dim=length(dim2))))
    dim1 = dim(array1)
  } else if(length(dim1) > length(dim2)){
    array2 = kronecker(array2, array(1, dim = array(1, dim=length(dim1))))
    dim2 = dim(array2)
  }

  if(min(dim1-dim2) < 0){
    array1_reshaped = kronecker(array1, array(1, dim=dim2-dim1+1))
    array2_reshaped = array2
  } else if(max(dim1-dim2) > 0){
    array1_reshaped = array1
    array2_reshaped = kronecker(array2, array(1, dim=dim1-dim2+1))
  } else{
    array1_reshaped = array1
    array2_reshaped = array2
  }
  if(operation == '*'){
    array_output = array1_reshaped*array2_reshaped
  } else if(operation == '/'){
    array_output = array1_reshaped/array2_reshaped
  } else if(operation == '+'){
    array_output = array1_reshaped+array2_reshaped
  } else if(operation == '-'){
    array_output = array1_reshaped-array2_reshaped
  } else if(operation == '^'){
    array_output = array1_reshaped^array2_reshaped
  }
  return(array_output)
}

#' Collapse array along one of the dimensions by adding the elements along that
#' dimension.
#'
#' @param array Array to collapse along one dimension.
#' @param dimension Dimension to collapse the array.
#'
#' @return An array that has been collapsed along the given dimension.
#'
#' @examples
sumDims = function(array, dimension){
  dim1 = length(dim(array))
  keep_dims = setdiff(seq(1, dim1, length.out = dim1), c(dimension))
  array_output = apply(array, MARGIN = keep_dims, FUN = 'sum')
  array_output = kronecker(array_output, array(1, dim=array(1, dim=dim1)))
  dim_reshape = seq(1, dim1, length.out = dim1)
  dim_reshape[(dimension+1):length(dim_reshape)] = dim_reshape[dimension:(length(dim_reshape)-1)]
  dim_reshape[dimension] = dim1
  array_output = aperm(array_output, dim_reshape);

  return(array_output)
}


#' Collapse array 2 along one of the dimensions by adding the elements along that
#' dimension.
#'
#' @param array Array to collapse along one dimension.
#' @param dimension Dimension to collapse the array.
#'
#' @return An array that has been collapsed along the given dimension.
sumDims2 = function(array, dimension){
  dim1 = length(dim(array))
  keep_dims = setdiff(seq(1, dim1, length.out = dim1), c(dimension))
  array_output = apply(array, MARGIN = keep_dims, FUN = 'sum')
  array_output = kronecker(array_output, array(1, dim=array(1, dim=dim1)))
  return(array_output)
}

