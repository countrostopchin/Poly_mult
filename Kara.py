def gf2_add(x, y):
  """Given two bit vectors x and y, returns their sum.
  Note: these need not be of the same size."""
  len_x = len(x)
  len_y = len(y)
  m = max(len_x, len_y)
  x = pad(x, m - len_x)
  y = pad(y, m - len_y)
  x_y = zip(x, y)
  return [ a ^ b for (a, b) in x_y ]  

def pad(x, i, bool=True):
  """Given a vector x of length n, an integer i and a Boolean
  bool (whose default value is True), either returns
  x padded on the left (resp., right) by i 0's
  if bool is True (resp., False)."""
  n = len(x)
  pad = [0]*i
  return pad + x if bool else x + pad

def gf2_const_prod(x, y):
  """Assuming that x and y are bit vectors, with at least one of them
  of length exactly 1, computes the product of x and y, viewed as big
  endian polynomials with binary coefficients."""
  len_x = len(x)
  len_y = len(y)
  min_poly = x if len_x < len_y else y
  max_poly = y if len_x < len_y else x
  min_len  = len(min_poly)
  max_len  = len(max_poly)
  if min_len==1:
    return max_poly if min_poly[0] else [0]*max_len
  else:
    pass

def gf2_poly_prod(x, y):
  """Given two bit vectors x and y, computes their product
  viewed as big endian polynomials with binary coefficients
  using the Karatsuba multiplication trick."""
  len_x = len(x)
  len_y = len(y)
  min_poly = x if len_x < len_y else y
  max_poly = y if len_x < len_y else x
  max_len  = len(max_poly)
  min_len = min(len_x, len_y)
  assert(min_len > 0)
  if min_len == 1:
    return gf2_const_prod(x, y)
  else:
    mid_len  = max_len // 2
    diff_len = max_len - mid_len
    max_p1  = max_poly[mid_len:]
    poly1   = gf2_poly_prod(max_p1, min_poly)
    max_p2  = max_poly[:mid_len]
    poly2   = gf2_poly_prod(max_p2, min_poly)
    poly2   = pad(poly2, diff_len, False)
    return gf2_add(poly1, poly2)