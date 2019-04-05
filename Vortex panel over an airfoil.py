"""
Created on Fri Jan 13 10:29:59 2017

author: Dipak Kumar Sisodiya (UR13AE021)

            ## Panel Method over an airfoil for calculation of lift ##
"""

# importing libraries and modules 
import numpy
from scipy import integrate, linalg
from matplotlib import pyplot
import math


# load geometry from data file using function Single_airfoil()
def Single_airfoil():
    print("input primary airfoil file location")

    file1 = input("Enter here: ")

    x1, y1 = numpy.loadtxt(file1, dtype=float, unpack=True)

    # using th angle of attack###
    xaoa = x1
    yaoa = y1
    return xaoa, yaoa

x,y = Single_airfoil()
AOB = float(input("Enter the Angle of attack:"))
con1 = (math.pi / 180)
AOA = AOB * con1
    
 
# plot geometry
val_x, val_y = 0.1, 0.2              
xp_min, xp_max = x.min(), x.max()    
yp_min, yp_max = y.min(), y.max()

# plot limits
xp_start, xp_end = xp_min-val_x*(xp_max-xp_min), xp_max+val_x*(xp_max-xp_min)
yp_start, yp_end = yp_min-val_y*(yp_max-yp_min), yp_max+val_y*(yp_max-yp_min)

# priliminary plotting the airfoil 
size = 10
pyplot.figure(figsize=(size, (yp_end-yp_start)/(xp_end-xp_start)*size))
pyplot.grid(True)
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(xp_start, xp_end)
pyplot.ylim(yp_start, yp_end)
pyplot.plot(x, y, color='k', linestyle='-', linewidth=2);

#### CLass named panle to compute panel number and size###################
class Panel:
    """Contains information related to a panel."""
    def __init__(self, xa, ya, xb, yb):
        """function to create a panel.
        
        Parameters
        ---------_
        xa, ya: float
            Coordinates of the starting-point.
        xb, yb: float
            Coordinates of the ending-point.
        """
        self.xa, self.ya = xa, ya                # panel starting-point
        self.xb, self.yb = xb, yb                # panel ending-point
        
        self.xc, self.yc = (xa+xb)/2, (ya+yb)/2         # panel center
        self.length = numpy.sqrt((xb-xa)**2+(yb-ya)**2) # panel length
        
        # orientation of panel (angle between x-axis and panel's normal)
        if xb-xa <= 0.0:
            self.beta = numpy.arccos((yb-ya)/self.length)
        elif xb-xa > 0.0:
            self.beta = numpy.pi + numpy.arccos(-(yb-ya)/self.length)
        
        # panel location
        if self.beta <= numpy.pi:
            self.loc = 'upper'                   # upper surface 
        else:
            self.loc = 'lower'                   # lower surface
        
        self.sigma = 0.0                         # source strength
        self.vt = 0.0                            # tangential velocity
        self.cp = 0.0                            # pressure coefficient

def define_panels(x, y, N=200):
    """Discretization of the geometry into panels using 'cosine' method.
    
    Parameters:
    x, y: Numpy 1d array (float)
        Coordinates of the geometry
    N: int
        Number of panels; default: 40.
    
    Returns:
    panels: Numpy 1d array (Panel object)
        Array of panels.
    """
    
    R = (x.max()-x.min())/2.0     # circle radius to be used for cosine method
    x_center = (x.max()+x.min())/2.0 # x-coordinate of circle center
    
    theta = numpy.linspace(0.0, 2.0*numpy.pi, N+1) # array of angles
    x_circle = x_center +  R*numpy.cos(theta)      # x-coordinates of circle
    
    x_ends = numpy.copy(x_circle)     # x-coordinate of panels end-points
    y_ends = numpy.empty_like(x_ends) # y-coordinate of panels end-points
    
    # extend coordinates to consider closed surface
    x, y = numpy.append(x, x[0]), numpy.append(y, y[0])
    
    # compute y-coordinate of end-points by projection y= ax + b
    I = 0
    for i in range(N):
        while I < len(x)-1:
            if (x[I] <= x_ends[i] <= x[I+1]) or (x[I+1] <= x_ends[i] <= x[I]):
            #Requires the value of x to be between the end values of x
                break
            else:
                I += 1
        a = ((y[I+1]-y[I])/(x[I+1]-x[I]))   # slope of the curve
        b = (y[I+1] - a*x[I+1])             # constant b value
        y_ends[i] = a*x_ends[i] + b       #  y = ax + b 
    y_ends[N] = y_ends[0]

    panels =numpy.empty(N ,dtype = object)
    for i in range (N):
       panels[i] = Panel(x_ends[i], y_ends[i], x_ends[i+1], y_ends[i+1])           
    return panels
    
# discretize geoemetry into panels
#N = int(input("Enter the number of panels : "))
panels = define_panels(x, y, N=50)

# plot discretized airofil geoemtry

val_x, val_y = 0.1, 0.2              

xp_min = min(panel.xa for panel in panels) #calling a object in class panel
xp_max = max(panel.xa for panel in panels)    
yp_min = min(panel.ya for panel in panels)
yp_max = max(panel.ya for panel in panels)    

# plot-limits
xp_start, xp_end = xp_min-val_x*(xp_max-xp_min), xp_max+val_x*(xp_max-xp_min)
yp_start, yp_end = yp_min-val_y*(yp_max-yp_min), yp_max+val_y*(yp_max-yp_min)

# plot
size = 10
pyplot.figure(figsize=(size, (yp_end-yp_start)/(xp_end-xp_start)*size))
pyplot.grid(True)
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(xp_start, xp_end)
pyplot.ylim(yp_start, yp_end)
pyplot.plot(x, y, color='k', linestyle='-', linewidth=2)
pyplot.plot(numpy.append([panel.xa for panel in panels], panels[0].xa), 
         numpy.append([panel.ya for panel in panels], panels[0].ya), 
         linestyle='-', linewidth=1, marker='o', markersize=6, color='#CD2305');

class Freestream:
    """Freestream conditions."""
    def __init__(self, u_inf=1.0, alpha=0.0):
        """Sets the freestream conditions.
        
        Parameters
        ----------
        u_inf: float
            Freestream speed; default: 1.0.
        alpha: float
            Angle of attack in degrees; default 0.0.
        """
        self.u_inf = u_inf
        self.alpha = alpha*numpy.pi/180.0 # degrees to radians

# define freestream conditions
freestream = Freestream(u_inf=20.0, alpha=AOB)

def integral(x, y, panel, dxdk, dydk):
    """Evaluates the contribution from a panel at a given point.
    
    Parameters
    ----------
    x, y: float
        Coordinates of the point.
    panel: Panel object
        Panel whose contribution is evaluated.
    dxdk: float
        Value of the derivative of x in a certain direction.
    dydk: float
        Value of the derivative of y in a certain direction.
    
    Returns
    -------
    Contribution from panel at point (x, y).
    """
    def func(s):
        return ( ((x - (panel.xa - numpy.sin(panel.beta)*s))*dxdk
                  +(y - (panel.ya + numpy.cos(panel.beta)*s))*dydk)
                / ((x - (panel.xa - numpy.sin(panel.beta)*s))**2
                   +(y - (panel.ya + numpy.cos(panel.beta)*s))**2) )
    return integrate.quad(lambda s:func(s), 0., panel.length)[0]

def source_contribution_normal(panels):
    """Builds the source contribution matrix for the normal velocity.
    
    Parameters
    ----------
    panels: Numpy 1d array (Panel object)
        List of panels.
    
    Returns
    -------
    A: Numpy 2d array (float)
        Source contribution matrix.
    """
    A = numpy.empty((panels.size, panels.size), dtype=float)
    # source contribution on a panel from itself
    numpy.fill_diagonal(A, 0.5)
    # source contribution on a panel from others
    for i, panel_i in enumerate(panels):
        for j, panel_j in enumerate(panels):
            if i != j:
                A[i, j] = 0.5/numpy.pi*integral(panel_i.xc, panel_i.yc, 
                                                panel_j,
                                                numpy.cos(panel_i.beta),
                                                numpy.sin(panel_i.beta))
    
    return A

def vortex_contribution_normal(panels):
    """Builds the vortex contribution matrix for the normal velocity.
    
    Parameters
    ----------
    panels: Numpy 1d array (Panel object)
        List of panels.
    
    Returns
    -------
    A: Numpy 2d array (float)
        Vortex contribution matrix.
    """
    A = numpy.empty((panels.size, panels.size), dtype=float)
    # vortex contribution on a panel from itself
    numpy.fill_diagonal(A, 0.0)
    # vortex contribution on a panel from others
    for i, panel_i in enumerate(panels):
        for j, panel_j in enumerate(panels):
            if i != j:
                A[i, j] = -0.5/numpy.pi*integral(panel_i.xc, panel_i.yc, 
                                                 panel_j,
                                                 numpy.sin(panel_i.beta),
                                                  -numpy.cos(panel_i.beta))
   
    return (A)

A_source = source_contribution_normal(panels)
B_vortex = vortex_contribution_normal(panels)

def kutta_condition(A_source, B_vortex):
    """Builds the Kutta condition array.
    
    Parameters
    ----------
    A_source: Numpy 2d array (float)
        Source contribution matrix for the normal velocity.
    B_vortex: Numpy 2d array (float)
        Vortex contribution matrix for the normal velocity.
    
    Returns
    -------
    b: Numpy 1d array (float)
        The left hand-side of the Kutta-condition equation.
    """
    b = numpy.empty(A_source.shape[0]+1, dtype=float)
    # matrix of source contribution on tangential velocity
    # is the same than
    # matrix of vortex contribution on normal velocity
    b[:-1] = B_vortex[0, :] + B_vortex[-1, :]
    # matrix of vortex contribution on tangential velocity
    # is the opposite of
    # matrix of source contribution on normal velocity
    b[-1] = - numpy.sum(A_source[0, :] + A_source[-1, :])
    return b
    
def build_singularity_matrix(A_source, B_vortex):
    """Builds the left hand-side matrix of the system
    arising from source and vortex contributions.
    
    Parameters
    ----------
    A_source: Numpy 2d array (float)
        Source contribution matrix for the normal velocity.
    B_vortex: Numpy 2d array (float)
        Vortex contribution matrix for the normal velocity.
    
    Returns
    -------
    A:  Numpy 2d array (float)
        Matrix of the linear system.
    """
    A = numpy.empty((A_source.shape[0]+1, A_source.shape[1]+1), dtype=float)
    # source contribution matrix
    A[:-1, :-1] = A_source
    # vortex contribution array
    A[:-1, -1] = numpy.sum(B_vortex, axis=1)
    # Kutta condition array
    A[-1, :] = kutta_condition(A_source, B_vortex)
    return A
    
def build_freestream_rhs(panels, freestream):
    """Builds the right hand-side of the system 
    arising from the freestream contribution.
    
    Parameters
    ----------
    panels: Numpy 1d array (Panel object)
        List of panels.
    freestream: Freestream object
        Freestream conditions.
    
    Returns
    -------
    b: Numpy 1d array (float)
        Freestream contribution on each panel and on the Kutta condition.
    """
    b = numpy.empty(panels.size+1,dtype=float)
    # freestream contribution on each panel
    for i, panel in enumerate(panels):
        b[i] = -freestream.u_inf * numpy.cos(freestream.alpha - panel.beta)
    # freestream contribution on the Kutta condition
    b[-1] = -freestream.u_inf*( numpy.sin(freestream.alpha-panels[0].beta)
                               +numpy.sin(freestream.alpha-panels[-1].beta) )
    return b
    
A = build_singularity_matrix(A_source, B_vortex)
b = build_freestream_rhs(panels, freestream)

# solve for singularity strengths
strengths = numpy.linalg.solve(A, b)

# store source strength on each panel
for i , panel in enumerate(panels):
    panel.sigma = strengths[i]
    
# store circulation density
gamma = strengths[-1]

def compute_tangential_velocity(panels, freestream, gamma, A_source, B_vortex):
    """Computes the tangential surface velocity.
    
    Parameters
    ----------
    panels: Numpy 1d array (Panel object)
        List of panels.
    freestream: Fresstream object
        Freestream conditions.
    gamma: float
        Circulation density.
    A_source: Numpy 2d array (float)
        Source contribution matrix for the normal velocity.
    B_vortex: Numpy 2d array (float)
        Vortex contribution matrix for the normal velocity.
    """
    A = numpy.empty((panels.size, panels.size+1), dtype=float)
    # matrix of source contribution on tangential velocity
    # is the same than
    # matrix of vortex contribution on normal velocity
    A[:, :-1] = B_vortex
    # matrix of vortex contribution on tangential velocity
    # is the opposite of
    # matrix of source contribution on normal velocity
    A[:, -1] = -numpy.sum(A_source, axis=1)
    # freestream contribution
    b = freestream.u_inf*numpy.sin([freestream.alpha-panel.beta 
                                    for panel in panels])
    
    strengths = numpy.append([panel.sigma for panel in panels], gamma)
    
    tangential_velocities = numpy.dot(A, strengths) + b
    
    for i, panel in enumerate(panels):
        panel.vt = tangential_velocities[i]

# tangential velocity at each panel center.
compute_tangential_velocity(panels, freestream, gamma, A_source, B_vortex)

def compute_pressure_coefficient(panels, freestream):
    """Computes the surface pressure coefficients.
    
    Parameters
    ----------
    panels: Numpy 1d array (Panel object)
        List of panels.
    freestream: Freestream object
        Freestream conditions.
    """
    for panel in panels:
        panel.cp = 1.0 - (panel.vt/freestream.u_inf)**2

# surface pressure coefficient
compute_pressure_coefficient(panels, freestream)

# plot surface pressure coefficient
val_x, val_y = 0.1, 0.2
x_min = min( panel.xa for panel in panels )
x_max = max( panel.xa for panel in panels )
cp_min = min( panel.cp for panel in panels )
cp_max = max( panel.cp for panel in panels )
x_start, x_end = x_min-val_x*(x_max-x_min), x_max+val_x*(x_max-x_min)
y_start, y_end = cp_min-val_y*(cp_max-cp_min), cp_max+val_y*(cp_max-cp_min)

pyplot.figure(figsize=(10, 6))
pyplot.grid(True)
pyplot.xlabel(r'$x$', fontsize=16)
pyplot.ylabel(r'$C_p$', fontsize=16)
pyplot.plot([panel.xc for panel in panels if panel.loc == 'upper'], 
         [panel.cp for panel in panels if panel.loc == 'upper'],
         label='upper surface',
         color='r', linestyle='-', linewidth=2, marker='o', markersize=6)
pyplot.plot([panel.xc for panel in panels if panel.loc == 'lower'], 
         [panel.cp for panel in panels if panel.loc == 'lower'], 
         label= 'lower surface',
         color='b', linestyle='-', linewidth=1, marker='o', markersize=6)
pyplot.legend(loc='best', prop={'size':14})
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.gca().invert_yaxis()
pyplot.title('Number of panels: {}'.format(panels.size), fontsize=16);
             
# calculate accuracy
accuracy = sum([panel.sigma*panel.length for panel in panels])
print('sum of singularity strengths: {:0.6f}'.format(accuracy))

# compute lift
cl = ( gamma*sum(panel.length for panel in panels)
       / (0.5*freestream.u_inf*(x_max-x_min)) )
print('lift coefficient: CL = {:0.3f}'.format(cl))

    
def get_velocity_field(panels, freestream, X, Y):
    """Returns the velocity field.
    
    Arguments
    ---------
    panels -- array of panels.
    freestream -- farfield conditions.
    X, Y -- mesh grid.
    """
    Nx, Ny = X.shape
    u, v = numpy.empty((Nx, Ny), dtype=float), numpy.empty((Nx, Ny), dtype=float)
    
    for i in range(Nx):
        for j in range(Ny):
            u[i,j] = freestream.u_inf*math.cos(freestream.alpha)\
                     + 0.5/math.pi*sum([p.sigma*integral(X[i,j], Y[i,j], p, 1, 0) for p in panels])
            v[i,j] = freestream.u_inf*math.sin(freestream.alpha)\
                     + 0.5/math.pi*sum([p.sigma*integral(X[i,j], Y[i,j], p, 0, 1) for p in panels])
    
    return u, v

# defines a mesh grid
Nx, Ny = 20, 20                  # number of points in the x and y directions
val_x, val_y = 1.0, 2.0
x_min, x_max = min( panel.xa for panel in panels ), max( panel.xa for panel in panels )
y_min, y_max = min( panel.ya for panel in panels ), max( panel.ya for panel in panels )
x_start, x_end = x_min-val_x*(x_max-x_min), x_max+val_x*(x_max-x_min)
y_start, y_end = y_min-val_y*(y_max-y_min), y_max+val_y*(y_max-y_min)

X, Y = numpy.meshgrid(numpy.linspace(x_start, x_end, Nx), numpy.linspace(y_start, y_end, Ny))

# computes the velicity field on the mesh grid
u, v = get_velocity_field(panels, freestream, X, Y)  
# plots the velocity field
size=10
pyplot.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.streamplot(X, Y, u, v, density=1, linewidth=1, arrowsize=1, arrowstyle='->')
pyplot.fill([panel.xc for panel in panels], 
         [panel.yc for panel in panels], 
         color='k', linestyle='solid', linewidth=2, zorder=2)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
# computes the pressure field
cp = 1.0 - (u**2+v**2)/freestream.u_inf**2

# plots the pressure field
size=12
pyplot.figure(figsize=(1.1*size, (y_end-y_start)/(x_end-x_start)*size))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
contf = pyplot.contourf(X, Y, cp, levels=numpy.linspace(-2.0, 1.0, 100), extend='both')
cbar = pyplot.colorbar(contf)
cbar.set_label('$C_p$', fontsize=16)
cbar.set_ticks([-2.0, -1.0, 0.0, 1.0])
pyplot.fill([panel.xc for panel in panels], 
         [panel.yc for panel in panels], 
         color='k', linestyle='solid', linewidth=2, zorder=2)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.title('Contour of pressure field');
