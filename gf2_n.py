import gf2

def gf2_add(p1, p2):
  """Assuming that p1 and p2 are vectors whose
  coefficients are bit vectors of length n, returns
  the sum of p1 and p2 (viewed as polynomials).
  Note: these need not be of the same size."""
  len_p1 = len(p1)
  len_p2 = len(p2)
  m = max(len_p1, len_p2)
  p1 = pad(x, m - len_p1)
  p2 = pad(y, m - len_p2)
  p1_p2 = zip(p1, p2)
  return [ gf2_add(a, b) for (a, b) in p1_p2 ]

def degree(f):
  """Assuming that f is a vector whose coefficients are
  bit vectors of length n, viewed as a big endian polynomial
  with coefficients in GF(2^n), returns the degree of f."""
  i = 0
  degree = len(f) - 1
  while f[i] == [0]*len(f[i]):
    i += 1
    degree -= 1
  return degree

def isZero(f):
  """Assuming that f is a vector whose coefficients are
  bit vectors of length n, viewed as a big endian polynomial
  with coefficients in GF(2^n), returns True if f is the
  zero polynomial and False else."""
  bool = degree(f) == -1
  return bool

def pad(p, i, n, bool=True):
  """Assuming that p is a vector whose coefficients are
  bit vectors of length n, viewed as a big endian polynomial
  with coefficients in GF(2^n), i is an integer and bool is a
  Boolean (whose default value is True), pads p on the left with
  i 0's if bool is True and on the right if bool is False."""
  pad = [[0] * n] * i
  return pad + p if bool else p + pad

def monomial(i, n):
  pass

def gf2n_const_prod(p1, p2, f):
   """Assuming that p1 and p2 are vectors whose
  coefficients are bit vectors of length n, viewed as big
  endian polynomials with coefficients in GF(2^n), that
  f is the modulus of GF(2^n), and that the length of at
  least one of p1 and p2 is exactly 1, returns product of
  p1 and p2."""
  len_p1 = len(p1)
  len_p2 = len(p2)
  min_poly = p1 if len_p1 < len_p2 else p2
  max_poly = p2 if len_p1 < len_p2 else p1
  min_len  = len(min_poly)
  max_len  = len(max_poly)
  if min_len==1:
    c = min_poly[0]
    p = [ gf2_poly_prod c y for y in max_poly ]
    return p
  else:
    pass

def gf2n_poly_prod(p1, p2, f):
   """Assuming that p1 and p2 are vectors whose
  coefficients are bit vectors of length n, viewed as big
  endian polynomials with coefficients in GF(2^n), and that
  f is the modulus of GF(2^n), returns the product of p1
  and p2."""
  len_p1 = len(p1)
  len_p2 = len(p2)
  min_poly = p1 if len_p1 < len_p2 else p2
  max_poly = p2 if len_p1 < len_p2 else p1
  min_len  = len(min_poly)
  max_len  = len(max_poly)
  assert(min_len > 0)
  if min_len == 1:
    return gf2n_const_prod(p1, p2, f)
  else:
    mid_len  = max_len // 2
    diff_len = max_len - mid_len
    max_p1  = max_poly[mid_len:]
    poly1   = gf2n_poly_prod(max_p1, min_poly)
    max_p2  = max_poly[:mid_len]
    poly2   = gf2n_poly_prod(max_p2, min_poly)
    poly2   = pad(poly2, diff_len, False)
    return gf2n_add(poly1, poly2)
