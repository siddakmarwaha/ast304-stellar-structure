########################################################################
# Team 6: Agrim, Patrick, and Siddak
# AST 304, Fall 2020
# Michigan State University
# This header (minus this line) should go at the top of all code files.
########################################################################

"""
<Description of this module goes here: what it does, how it's used.>
"""

# all routines that take a single step should have the same interface
# fEuler is complete, except for documentation. you can use this as a pattern 
# for the other two routines.
def fEuler(f,t,z,h,args=()):
    """
    <Description of routine goes here: what it does, how it's called>
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        <fill this in>
    
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """
    
    # The following trick allows us to pass additional parameters to f
    # first we make sure that args is of type tuple; if not, we make it into
    # that form
    if not isinstance(args,tuple):
        args = (args,)
    
    # when we call f, we use *args to pass it as a list of parameters.
    # for example, if elsewhere we define f like
    # def f(t,z,x,y):
    #    ...
    # then we would call this routine as
    # znew = fEuler(f,t,z,h,args=(x,y))
    #
    return z + h*f(t,z,*args)

# You will need to flesh out the following routines for a second-order
# Runge-Kutta step and a fourth order Runge-Kutta step.

def rk2(f,t,z,h,args=()):
    """
    <Description of routine goes here: what it does, how it's called>
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        <fill this in>
    
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """
    
    if not isinstance(args,tuple):
        args = (args,)
    
    # delete the line "pass" when you put in the full routine
    # Step 1: Calculate the midpoint using the forward Euler step
    zP = z + (h/2)*f(t,z,*args)
    
    # Step 2: Use the midpoint estimate to compute the final step
    znew = z + h*f(t+h/2,zP,*args)
    return znew

def rk4(f,t,z,h,args=()):
    """
    <Description of routine goes here: what it does, how it's called>
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        <fill this in>
    
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """
   
    if not isinstance(args,tuple):
        args = (args,)
    
    # delete the line "pass" when you put in the full routine
    # Step 1: Calculate k1
    k1 = f(t,z,*args)
    
    # Step 2: Calculate k2 using the midpoint t + h/2 and the slope k1
    k2 = f(t+h/2,z+(h/2)*k1,*args)
    
    # Step 3: Calculate k3, another midpoint estimate using k2
    k3 = f(t+h/2,z+(h/2)*k2,*args)
    
    # Step 4: Calculate k4 at the endpoint t + h using k3
    k4 = f(t+h,z+h*k3,*args)
    
    # Final calculation: weighted sum of k1, k2, k3, k4
    znew = z + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    return znew
