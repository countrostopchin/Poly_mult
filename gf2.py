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

def degree(f):
  """Given a bit vector f, viewed as a big endian polynomial
  with binary coefficients, returns the degree of f."""
  i = 0
  degree = len(f) - 1
  while not f[i]:
    i += 1
    degree -= 1
  return degree

def isZero(f):
  """Given a bit vector f, viewed as a big endian polynomial
  with binary coefficients, returns True if f is identically
  zero and True else."""
  bool = degree(f)== -1
  return bool

def pad(x, i, bool=True):
  """Given a vector x of length n, an integer i and a Boolean
  bool (whose default value is True), either returns
  x padded on the left (resp., right) by i 0's
  if bool is True (resp., False)."""
  n = len(x)
  pad = [0]*i
  return pad + x if bool else x + pad

def monomial(i, n):
  """Given integers i and n, returns a bitvector of length n
  all of whose entries are False, with the exception of n - i."""
  monom = [0]*n
  assert(i>0)
  monom[n-i] = 1
  return monom

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

def gf2_rem(x, f):
  """Given bitvectors x and f, viewed as big
  endian polynomials with binary coefficients,
  returns the result of dividing (with remainder)
  the polynomial x by f. The output is the
  remainder."""
  assert(not isZero(f))
  len_x = len(x)
  len_f = len(f)
  deg_f = degree(f)
  if len_x < deg_f:
    return x
  rem = x
  while degree(rem) >= degree(f):
    diff  = degree(rem) - degree(f)
    monom = monomial(diff, len_x)
    prod  = gf2_prod(monom, f)   ## need to replace this with a left shift operator
    rem   = gf2_add(rem, prod)
  return rem[-len_f+1:]

def gf2_div(x, f):
  """Given bitvectors x and f, viewed as big
  endian polynomials with binary coefficients,
  returns the result of dividing (with remainder)
  the polynomial x by f. The output is the
  quotient."""
  assert(not isZero(f))
  len_x = len(x)
  len_f = len(f)
  deg_f = degree(f)
  if len_x < deg_f:
    return x
  quot = [0]*len_x
  rem  = x
  while degree(rem) >= degree(f):
    diff  = degree(rem) - degree(f)
    monom = monomial(diff, len_x)
    quot  = gf2_add(quot, monom)
    prod  = gf2_prod(monom, f)   ## need to replace this with a left shift operator
    rem   = gf2_add(rem, prod)
  return quot

def gf2_mod_prod(x,y,f):
  """Given three bit vectors x, y and f, viewed as big
  endian polynomials with binary coefficients, computes
  the remainder of x * y when divided by f."""
  return gf2_rem(gf2_poly_prod(x, y), f)

def gf2_gcd(x, y):
  """Given two bit vectors x and y, viewed as big endian polynomials
  with binary coefficients, computes the gcd of x and y."""
  if isZero(y):
    return x
  else:
    r = gf2_rem(x, y)
    return gf2_gcd(y, r)

def gf2_rel_prime(x, y):
  """Given two bit vectors x and y, viewed as big endian polynomials
  with binary coefficients, returns True if x and y are relatively
  prime and False else."""
  bool = degree(gf2_gcd(x, y)) == 0
  return bool

def gf2_inv(x, y):
  """Given two bit vectors x and y, viewed as big endian polynomials
  with binary coefficients, computes the inverse of x modulo y,
  assuming that x and y are relatively prime."""
  assert(gf2_rel_prime(x, y))
  len_y = len(y)
  b = [0]*len_y
  a = [0]*(len_y - 1) + [1]
  x = x
  y = y
  r = gf2_rem(x, y)
  while degree(r) >= 0:
    b1   = b
    quot = gf2_div(x, y)
    prod = gf2_poly_prod(quot, b1)
    b    = gf2_add(a, prod)
    a    = b1
    x    = y
    y    = r
    r    = gf2_rem(x, y)
  return b




  
