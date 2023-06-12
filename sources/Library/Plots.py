from matplotlib.pyplot import figure, plot, bar, scatter, subplot, title, xlabel, ylabel, text, axis, show, legend, grid, contour, contourf, axes
from numpy import polyfit, poly1d, zeros, shape, reshape, ceil, absolute, amax, transpose
from matplotlib.animation import FuncAnimation


""" Principal Functions """
def Figure():
    return figure()

def Plot(X,Y,Title=[],Xlabel=[],Ylabel=[],Text=[],Axis=[],Show=[],Form=[],**arg):
    
    p = plot(X,Y, Form, **arg)
    if Title: title(Title)
    if Xlabel: xlabel(Xlabel)
    if Ylabel: ylabel(Ylabel)
    if Text: text(Text[0],Text[1],Text[2])
    if Axis: axis([Axis[0],Axis[1],Axis[2],Axis[3]])
    if not Show or Show=="show": show()
    
    return p

def Fig3D(**args):
    ax = axes(projection='3d')
    return ax


# def subplots(XXX,**kwargs): # XXX = [nrows,ncols,index] or nrowsncolsindex
    
#     # In Progress
#     subplot(XXX[0], XXX[1], XXX[2], **kwargs) if len(XXX) == 3 else subplot(XXX, **kwargs)
    
#     return

def Trendline(X,Y,**arg):
    
    p = poly1d(polyfit(X,Y, 2))
    
    return Plot(X,p(X),**arg)

def Bar(X,Y,Title=[],Xlabel=[],Ylabel=[],Text=[],Axis=[],Show=[]):
    
    b = bar(X,Y)
    if Title: title(Title)
    if Xlabel: xlabel(Xlabel)
    if Ylabel: ylabel(Ylabel)
    if Text: text(Text[0],Text[1],Text[2])
    if Axis: axis([Axis[0],Axis[1],Axis[2],Axis[3]])
    if not Show or Show=="show": show()

    return b
    
def Scatter(X,Y,Title=[],Xlabel=[],Ylabel=[],Text=[],Axis=[],Show=[],Form=[],**arg):
    
    s = scatter(X,Y, Form, **arg)
    if Title: title(Title)
    if Xlabel: xlabel(Xlabel)
    if Ylabel: ylabel(Ylabel)
    if Text: text(Text[0],Text[1],Text[2])
    if Axis: axis([Axis[0],Axis[1],Axis[2],Axis[3]])
    if not Show or Show=="show": show()
    
    return s

def Scatter3D(X,Y,Z,**args):
    ax = axes(projection='3d')
    ax.scatter3D(X,Y,Z,**args)
    return ax

def Contour(X,Y,Z, ctype = "lines", Show = [], **arg):

    if not ctype or ctype == "lines":
        c = contour(X, Y, Z, **arg)
    elif ctype == "filled":
        c = contourf(X, Y, Z, **arg)
    if not Show or Show=="show": show()

    return c

def Contour3D(X,Y,Z, ctype = "lines",**args):
    ax = axes(projection='3d')
    if not ctype or ctype == "lines":
        ax.contour(X,Y,Z,**args)
    elif ctype == "filled":
        ax.contourf(X,Y,Z,**args)
    return ax
    
def update_plot(num, walks, lines):
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[:num, :2].T)
        line.set_3d_properties(walk[:num, 2])
    return lines

def paint(N,Nc,r,index):
  line=zeros((N,Nc))
  for i in range(N):
    line[i,:]=r[i,index,:]
  return line

def animate_plot(r):

    N,Nb,Nc = shape(r)
    
    # Data: Nb lines as (num_steps, 3) arrays
    datas = [paint(N,Nc,r,index) for index in range(Nb)]

    # Attaching 3D axis to the figure
    fig = figure()
    ax = fig.add_subplot(projection="3d")

    # Create lines initially without data
    lines = [ax.plot([], [], [])[0] for _ in datas]

    # Setting the axes properties
    ax.set(xlim3d=(-10, 10), xlabel='X')
    ax.set(ylim3d=(-10, 10), ylabel='Y')
    ax.set(zlim3d=(-10, 10), zlabel='Z')

    # Creating the Animation object
    ani = FuncAnimation(fig, update_plot, N//2, fargs=(datas, lines), interval=2500)

    show()
    
    return ani


# --------------------------------------------------------------------------------- #


""" Complementary Functions """

def Xlabel(lab):
    return xlabel(lab)

def Ylabel(lab):
    return ylabel(lab)

def Legend(**arg):
    
    return legend(**arg)

def ShowPlot():
    
    return show()

def Grid():
    
    return grid()

def ClearPlot():
    
    #In Progress
    
    return

    