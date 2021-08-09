def title_page(plt,titles,name=None,paper=None):
    import matplotlib.pyplot as plt   
    pheight = 11  # page height, inches
    pwidth = 8.5  # page weight, inches

    fig = plt.figure(figsize=(pwidth,pheight))    
    
    ytop = 1.0/pheight  # Top margin, 1 inch, normalized
    xleft = 2.0/pwidth  # Left margin, 1 inch, normalized
    xright=1-2.5/pwidth
    textheight = 0
    spacing = 0

    main_title='Rayleigh Simulation Summary\n\n'
    if (paper != None):
        main_title+=paper+'\n\n'
    
    if (name != None):
        main_title += 'Simulation: '+name+'\n'
    
    fig.text(0.5,1-ytop,main_title,ha='center',va='top')
    
    clen = 72  # Number of characters per line in title page
    

    yline = 1-ytop*3

    left = 'Figure\n\n'
    right = 'Page\n\n'
    
    for i,t in enumerate(titles):
        left+=t+'\n'
        pnum = str(i+2)
        right+=pnum+'\n'
            
    fig.text(xleft,yline,left,va='top')
    fig.text(xright,yline,right,va='top')
