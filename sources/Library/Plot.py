from matplotlib.pyplot import Figure, plot, bar, title, xlabel, ylabel, text, axis, show



def Plot(X,Y,Title=[],Xlabel=[],Ylabel=[],Text=[],Axis=[],Show=[]):
    
    Figure
    plot(X,Y)
    if Title: title(Title)
    if Xlabel: xlabel(Xlabel)
    if Ylabel: ylabel(Ylabel)
    if Text: text(Text[0],Text[1],Text[2])
    if Axis: axis([Axis[0],Axis[1],Axis[2],Axis[3]])
    if not Show or Show=="show": show()
    
    return

def Bar(X,Y,Title=[],Xlabel=[],Ylabel=[],Text=[],Axis=[],Show=[]):
    
    Figure
    bar(X,Y)
    if Title: title(Title)
    if Xlabel: xlabel(Xlabel)
    if Ylabel: ylabel(Ylabel)
    if Text: text(Text[0],Text[1],Text[2])
    if Axis: axis([Axis[0],Axis[1],Axis[2],Axis[3]])
    if not Show or Show=="show": show()
    
    return

def clearplot():
    
    return