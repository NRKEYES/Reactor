def shotgun_graph(x,y, xlabel = '', ylabel= ''):
    fig = plt.figure()
    fig.set_size_inches(20, 8.5)
    title = "Shotgun Error Plot (kJ/mol)"
    fig.suptitle(title, fontsize = 25)
    plt.ylabel(ylabel, fontsize = 30)
    plt.xlabel(xlabel, fontsize = 30)
    
    plt.axhline(0,ls = "--", color = sb.xkcd_rgb["dusty purple"])
    plt.axvline(0,ls = "--", color = sb.xkcd_rgb["dusty purple"])
    
    ax = fig.add_subplot(111)
    ax.yaxis.label.set_color(sb.xkcd_rgb["dusty purple"])
    ax.xaxis.label.set_color(sb.xkcd_rgb["dusty purple"])
    ax.tick_params(axis='y', colors=sb.xkcd_rgb["dusty purple"], labelsize = 20)
    ax.tick_params(axis='x', colors=sb.xkcd_rgb["dusty purple"], labelsize = 20)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    
    plt.scatter(x,y, s = 100)
    plt.show()
    fig.savefig(str(xlabel+ylabel+'.png'))


    def stateLines(index, xSpacing = 5, stateSize = 2):
        # Blank lists
    verts = []
    drawState = []
    
    for key, state in myReaction.reaction_states.items():
        # Get vertices x and y
        verts.append(((np.floor(float(state.key)/1)*xSpacing),
                         kJ.get_value(index,state.key)))

        verts.append(((np.floor(float(state.key)/1)*xSpacing)+stateSize,
                         kJ.get_value(index,state.key)))
    
        # Move the pen and draw a line between the vertices we just made
        drawState.append(Path.MOVETO)
        drawState.append(Path.LINETO)
        
    return verts, drawState




def connectLines(index, xSpacing = 5, stateSize = 2 ):
    verts = []
    instructions = []
    
    
    for key, state in myReaction.reaction_states.items():
        if len(state.goesTo) > 0:
            for nextState in state.goesTo:
                verts.append(((np.floor(float(state.key)/1)*xSpacing)+stateSize, #x
                             kJ.get_value(index,state.key)))        # Y
                verts.append(((np.floor(float(state.key)/1+1)*xSpacing),
                             kJ.get_value(index,myReaction.reaction_states[str(nextState)].key)))  
                instructions.append(Path.MOVETO)
                instructions.append(Path.LINETO)
                     
    return verts, instructions




def subway_graph(index):

    xSpacing = 6
    stateSize = 2

    stateVert, stateRules = stateLines(index, xSpacing=xSpacing)
    path = Path(stateVert, stateRules)
    patch = patches.PathPatch(path,color = '#00CED1', lineWidth=7, capstyle = 'round')

    connectVert, connectRules = connectLines(index, xSpacing=xSpacing)
    path2 = Path(connectVert, connectRules)
    patch2 = patches.PathPatch(path2, color =  'lavender', lineWidth=7, capstyle = 'round')




    fig = plt.figure()
    fig.set_size_inches(20, 8.5)
    title = str(index[0]) + " Functional with " + index[1] + " Basis Set"
    fig.suptitle(title, fontsize = 30)
    plt.ylabel("kJ/mol", fontsize = 30)

    ax = fig.add_subplot(111)
    ax.yaxis.label.set_color(sb.xkcd_rgb["dusty purple"])
    ax.tick_params(axis='y', colors=sb.xkcd_rgb["dusty purple"],labelsize = 20)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().set_visible(False)
    
    ax.add_patch(patch2)
    ax.add_patch(patch)
    
    plt.axhline(0,ls = "--", color = sb.xkcd_rgb["dusty purple"])

    
    for key, state in myReaction.reaction_states.items():
        ax.text(float(state.key)*xSpacing, kJ.get_value(index,state.key)-18 , state.key,
                #verticalalignment='bottom', horizontalalignment='right',
                #transform=ax.transAxes,
                color=sb.xkcd_rgb["denim blue"], fontsize=10)

    xs, ys = zip(*stateVert)
    ax.plot(xs, ys, "o",lw=0, color='black', ms=7)
    ax.plot(xs, ys, "o",lw=0, color='white', ms=5)
   

    plt.show()
    fig.savefig(str(index[0])+str(index[1]+'.png'))
    
    
    
for index, row in kJ.iterrows():
    
    try:
        subway_graph(index)
    except:
        pass
    #break # to limit to one graph while troubleshooting
