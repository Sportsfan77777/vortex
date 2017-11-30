

###############################################################################

def load_data(beam):
    """ load contrasts for a particular beam size """
    directory = "beam%03d/%s" % (beam, save_directory)

    save_fn = "%s/id%04d_contrasts_lambda%04d_beam%03d.p" % (save_directory, id_number, wavelength, beam)
    contrasts = pickle.load(open(save_fn, "wb"))

    return contrasts


###############################################################################

##### PLOTTING #####

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

def make_plot(show = False):
    # Set up figure
    fig = plot.figure(figsize = (7, 6), dpi = dpi)
    ax = fig.add_subplot(111)

    ### Line Plots ###

    for i, beam in enumerate(beam_sizes):
        contrasts = load_data(beam)

        x = frame_range
        y = np.array(contrasts)
        plot.plot(x, y, c = colors[i], linewidth = linewidth, label = r"$%d$ $\mathrm{AU}$" % beam)

    # Annotate
    plot.xlabel("Time (planet orbits)", fontsize = fontsize)
    plot.ylabel("Contrasts", fontsize = fontsize)
    #plot.title("")

    plot.legend(loc = "upper left")

    # Axes
    plot.xlim(frame_range[0], frame_range[-1])
    if max_y is not None:
        plot.ylim(0, max_y)

    # Save, Show, and Close
    if version is None:
        save_fn = "%s/id%04d_comparing_contrasts_lambda%04d.png" % (save_directory, id_number, wavelength)
    else:
        save_fn = "%s/v%04d_id%04d_comparing_contrasts_lambda%04d_beam%03d.png" % (save_directory, version, id_number, wavelength)
    plot.savefig(save_fn, bbox_inches = 'tight', dpi = dpi)

    if show:
        plot.show()

    plot.close(fig) # Close Figure (to avoid too many figures)


###############################################################################

##### Make Plots! #####

make_plot(show = show)
