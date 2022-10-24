from matplotlib.pyplot import Figure, plot, bar, scatter, subplot, title, xlabel, ylabel, text, axis, show, legend
from numpy import polyfit, poly1d

""" Principal Functions """

def Plot(X,Y,Title=[],Xlabel=[],Ylabel=[],Text=[],Axis=[],Show=[],Form=[],**arg):
    
    Figure
    p = plot(X,Y, Form, **arg)
    if Title: title(Title)
    if Xlabel: xlabel(Xlabel)
    if Ylabel: ylabel(Ylabel)
    if Text: text(Text[0],Text[1],Text[2])
    if Axis: axis([Axis[0],Axis[1],Axis[2],Axis[3]])
    if not Show or Show=="show": show()
    
    return p

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
    
    Figure
    s = scatter(X,Y, Form, **arg)
    if Title: title(Title)
    if Xlabel: xlabel(Xlabel)
    if Ylabel: ylabel(Ylabel)
    if Text: text(Text[0],Text[1],Text[2])
    if Axis: axis([Axis[0],Axis[1],Axis[2],Axis[3]])
    if not Show or Show=="show": show()
    
    return s


""" Complementary Functions """

def Legend(**arg):
    
    return legend(**arg)

def ShowPlot():
    
    return show()

def ClearPlot():
    
    #In Progress
    
    return