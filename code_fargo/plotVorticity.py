"""
plot 2-D density maps

python plotVorticity.py
python plotVorticity.py frame_number
python plotVorticity.py -m
"""

save_directory = "vorticityMaps"

### Movie Commands ###
def make_movies():
    # Movie Parameters
    fps = 5

    path = save_directory + "/vorticityMap_%03d.png"
    output = save_directory + "/vorticityMap.mov"

    zoom_path = save_directory + "/zoom_vorticityMap_%03d.png"
    zoom_output = save_directory + "/vorticityMap_zoom.mov"

    # Movie Command
    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
    split_command = command.split()
    subprocess.Popen(split_command)

    command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, zoom_path, zoom_output)
    split_command = command.split()
    subprocess.Popen(split_command)

# Make only movies and then return
if (len(sys.argv) > 1) and (sys.argv[1] == "-m"):
    make_movies()
    # Terminate
    quit()


### Get FARGO Parameters ###
# Create param file if it doesn't already exist
param_fn = "params.p"
if not os.path.exists(param_fn):
    command = "python pickleParameters.py"
    split_command = command.split()
    subprocess.Popen(split_command)
fargo_par = pickle.load(open(param_fn, "rb"))

num_rad = np.loadtxt("dims.dat")[-2]
num_theta = np.loadtxt("dims.dat")[-1]

rad = np.loadtxt("used_rad.dat")[:-1]
theta = np.linspace(0, 2 * np.pi, num_theta)

surface_density_zero = float(fargo_par["Sigma0"])
scale_height = float(fargo_par["AspectRatio"])



##### PLOTTING #####

# Make Directory
try:
    os.mkdir(save_directory)
except:
    print "Directory Already Exists"

# Plot Parameters
cmap = "RdYlBu_r"
clim = [0, 2]

fontsize = 14
my_dpi = 100



##### Plot One File or All Files #####

if len(sys.argv) > 1:
    frame_number = int(sys.argv[1])
    make_plot(frame_number)
else:
    # Search for maximum frame
    density_files = glob.glob("gasvrad*.dat")
    max_frame = 0
    for d_f in density_files:
        name = d_f.split(".")[0] # for "gasdens999.dat", just "gasdens999"
        frame_number = int(name[7:]) # just 999
        if frame_number > max_frame:
            max_frame = frame_number
    num_frames = max_frame + 1

    #for i in range(num_frames):
    #    make_plot(i)

    p = Pool(8)
    p.map(make_plot, range(num_frames))
    p.terminate()

    #### Make Movies ####
    make_movies()


