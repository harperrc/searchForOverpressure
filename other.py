import math
import numpy as np

def looseLabel(min_val,max_val,**kwargs):

#  p.heckbert (12/02/1998) Graphics Gems, Academic Press 1990
#  modifications by c.harper

#  give min,max and optional # of steps will come up with
#  pretty label tics and markers

   ntick  = 5
   fmtStr = '%d'

   retVal = ypstruct()

   for name,value in kwargs.items():
      unknown = 1
      name_l  = name.lower()

      if (name_l == "ntick"):
         unknown = 0
         ntick   = value;

      if (name_l == "fmtstr"):
         unknown = 0
         fmtStr  = value;

      if (unknown):
         print('unknown key %s' %(name))

   rng      = nicenum(max_val - min_val,0)
   d        = nicenum(rng / (ntick - 1),1)
   graphmin = math.floor(min_val / d) * d
   graphmax = math.ceil(max_val / d) * d

   retVal.ntick   = int((graphmax - graphmin) / d) + 1;

   retVal.min_val = graphmin;
   retVal.max_val = graphmax;

   retVal.stp_val = d;

   retVal.tick_values = np.linspace(graphmin,graphmax,retVal['ntick']);
   retVal.tick_labels = []

   val = graphmin;
   for k in range(0,retVal.ntick):
      retVal.tick_labels.append(fmtStr % (val))
      val += d

   return (retVal)

def nicenum(x,round_dir):

#  old fortran code to round number to 'nice' values

   expv = math.floor(math.log10(x))

   f    = x / (10.0**expv)

   if (round_dir):
      if (f < 1.5):
         nf = 1.0
      elif (f < 3.0):
         nf = 2.0
      elif (f < 7.0):
         nf = 5.0
      else:
         nf = 10
   else:
      if (f < 1.0):
         nf = 1.0
      elif (f <= 2.0):
         nf = 2.0
      elif (f <= 5.0):
         nf = 5.0
      else:
         nf = 10

   value = nf * 10.0**expv

   return value

class ypstruct(dict):

# retreived from https://github.com/smkalami/ypstruct
# changed name struct to ypstruct to avoid any confusion (c. harper)
# Import copy as _copy, to perform deep copy operations (moved to main header)

    """
    A class to implement the C/C++ or MATLAB-like structures
    """
    
    def __repr__(self):
        """
        String representation of the struct
        """
        return "ypstruct({})".format(super().__repr__())

    def __getattr__(self, field):
        """
        Gets value of a field
        """
        if field not in dir(self):
            if field in self.keys():
                return self[field]
            else:
                return None
        else:
            return None
    
    def __setattr__(self, field, value):
        """
        Sets value of a field
        """
        if field not in dir(self):
            self[field] = value
        else:
            return super().__setattr__(field, value)
    
    def fields(self):
        """
        Gets the list of defined fields of the struct
        """
        return list(self.keys())

    def remove_field(self, field):
        """
        Removes a field from the struct
        """
        if field in self.keys():
            del self[field]
    
    def add_field(self, field, value = None):
        """
        Adds a new field to the struct
        """
        if field not in self.keys():
            self[field] = value

    def copy(self):
        """
        Creates a shallow copy of the struct
        """
        self_copy = ypstruct()
        for field in self.keys():
            if isinstance(self[field], ypstruct):
                self_copy[field] = self[field].copy()
            else:
                self_copy[field] = _copy.copy(self[field])
        
        return self_copy

    def deepcopy(self):
        """
        Creates a deep copy of the struct
        """
        self_copy = ypstruct()
        for field in self.keys():
            if isinstance(self[field], ypstruct):
                self_copy[field] = self[field].deepcopy()
            else:
                self_copy[field] = _copy.deepcopy(self[field])
        
        return self_copy

    def repeat(self, n):
        """
        Repeats/replicates the struct to create an array of structs (eg. for initialization)
        """
        return [self.deepcopy() for i in range(n)]

    def __mul__(self, n):
        """
        Overload * operator (multiplication) to repeat/replicate the struct
        """
        if not isinstance(n, int) and not isinstance(n, float):
            raise TypeError("Only integers are allowed.")
        return self.repeat(n)

    def __add__(self, other):
        """
        Overload + operator (addition) to merge two struct objects
        """
        if not isinstance(other, dict):
            raise TypeError("Only structure and dict objects are allowed.")
        result = self.deepcopy()
        result.update(other)
        return result

    def print(self):

       keys = self.fields()
       nc   = [len(i) for i in keys]
       mxnc = max(nc) + 1

       nv = 0

       for aKey in keys:
          print('%s' % (aKey),end='')
          for i in range(0,mxnc-nc[nv]):
             print(' ',end='')

          theType = type(self[aKey])

          if (theType == int):
             print('%20d' % (self[aKey]))
          elif (theType == bool):
             print('%r' % (self[aKey]))
          elif (theType == float or theType == np.float64):
             if (abs(self[aKey]) <= 1.0e-6):
                print('%20.8e' % (self[aKey]))
             elif (abs(self[aKey]) >= 1.0e6):
                print('%20.8e' % (self[aKey]))
             else:
                print('%20.12f' % (self[aKey]))
          elif (theType == str):
             print('%s' % (self[aKey]))
          elif (theType == list):
             nx = len(self[aKey])
             theType = type(self[aKey][0])
             if (theType == int):
                for j in range(0,nx):
                   print('%d ' % (self[aKey][j]),end='')
                print()
             elif (theType == bool):
                for j in range(0,nx):
                   print('%r ' % (self[aKey][j]),end='')
                print()
             elif (theType == float or theType == np.float64):
                for j in range(0,nx):
                   print('%f ' % (self[aKey][j]),end='')
                print()

          nv = nv + 1
