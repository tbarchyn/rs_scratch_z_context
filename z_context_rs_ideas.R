# Point cloud based remote sensing scratch
# copyright 2016 tom barchyn
# this has no warranty whatsoever!

# purpose: scratch sheet for working with point based remote sensing

# Notes:
# from visualsfm:
# save point cloud with colors as ply file (save NView match -> change 
#    extension to ply)
# open ply file in meshlab
# crop / rotate / scale, etc. (freeze rotation sets the coordinates)
# save as a JSON file format (yes, this is not a good idea, but it
#     works for now!)

library (rjson)     #jsonlite is too heavy duty for this
library (raster)
library (rgdal)

setwd ("C://Users//tom//Desktop//legoscene//models")

####################################################################
# read in the file and pull out raw vertices and associated colors
d <- fromJSON (file = "dense_cut.json")
vts <- d$vertices[[1]][5]
vts <- vts[[1]]
cols <- d$vertices[[3]][5]
cols <- cols[[1]]

# ok, now we can reassemble them
n <- length(vts) / 3
verts <- data.frame (id = 1:n)
verts$x <- vts[rep(c(T,F,F), n)]
verts$y <- vts[rep(c(F,T,F), n)]
verts$z <- vts[rep(c(F,F,T), n)]
verts$r <- cols[rep(c(T,F,F,F), n)]
verts$g <- cols[rep(c(F,T,F,F), n)]
verts$b <- cols[rep(c(F,F,T,F), n)]
verts$a <- cols[rep(c(F,F,F,T), n)]

# write it out
write.csv (verts, 'pointcloud.csv')


####################################################################
# CALC METRIC FROM THE POINT CLOUD
# first make that verts dataframe a spatial points dataframe
verts <- SpatialPointsDataFrame(data = verts, coords = cbind(verts$x, verts$y))

# make an overlay raster for aggregation purposes (raster object can mirror on disk)
csize <- 0.02
x_crds <- seq (min(verts$x), max(verts$x), by = csize)
y_crds <- seq (min(verts$y), max(verts$y), by = csize)
init_matrix <- matrix(NA, nrow = (length(y_crds)), ncol = length(x_crds))
protoras <- raster (init_matrix,
                    xmn = x_crds[1] - csize/2, xmx = x_crds[length(x_crds)] + csize/2,
                    ymn = y_crds[1] - csize/2, ymx = y_crds[length(y_crds)] + csize/2)

subset_verts <- function (verts, x_crds, y_crds, csize, i, j) {
    # function to subset the points from the verts dataframe
    # note: this won't scale, need to use postgis or something faster here for the
    # query, need something that doesn't need to hold everything in memory
    # returns subset of verts that are within the defined cell
    # i = the row
    # j = the col
    
    minx <- x_crds[j] - csize / 2
    maxx <- x_crds[j] + csize / 2
    miny <- y_crds[length(y_crds) + 1 - i] - csize / 2
    maxy <- y_crds[length(y_crds) + 1 - i] + csize / 2
    
    sub_verts <- verts [(verts$x > minx & verts$x <= maxx & verts$y > miny &
                             verts$y <= maxy), ]
    return (sub_verts)
}

# now make some rasters to store our outputs
mean_z <- protoras          # mean z val from all verts
max_z <- protoras
min_z <- protoras
mean_r <- protoras          # mean r from all verts
mean_g <- protoras
mean_b <- protoras
numverts <- protoras        # the number of verts
top_r <- protoras           # the r from the top point in the subset
top_g <- protoras
top_b <- protoras
bot_r <- protoras           # the r from the bottom point in the subset
bot_g <- protoras
bot_b <- protoras
mean_r_upper <- protoras    # metrics from the upper half of the points
mean_r_lower <- protoras    # metrics from the lower half of the points
mean_g_upper <- protoras
mean_g_lower <- protoras
mean_b_upper <- protoras
mean_b_lower <- protoras

# Loop over rasters
for (i in 1:nrow(protoras)) {
    for (j in 1:ncol(protoras)) {
        # subset the verts
        sverts <- subset_verts (verts, x_crds, y_crds, csize, i, j)
        
        # only calc metrics if we actually pulled some points
        if (nrow (sverts) > 0) {
            # calc the whole sample metrics (note this will crash with any NAs!)
            mean_z[i, j] <- mean (sverts$z)
            max_z[i, j] <- max (sverts$z)
            min_z[i, j] <- min (sverts$z)
            mean_r[i, j] <- mean (sverts$r)
            mean_g[i, j] <- mean (sverts$g)
            mean_b[i, j] <- mean (sverts$b)
            numverts[i, j] <- nrow (sverts)
            
            # get the top point and associated metrics (only if we get 1 return,
            # improvement would be handling weird cases with 2+ returns for min equality,
            # or just not do this if there is only 1 point in sverts, . . .)
            sverts_top <- sverts [sverts$z == max(sverts$z), ]
            if (nrow(sverts_top) == 1) {
                top_r[i, j] <- sverts_top$r[1]
                top_g[i, j] <- sverts_top$g[1]
                top_b[i, j] <- sverts_top$b[1]
            }
            
            # get the bottom point and associated metrics
            sverts_bot <- sverts [sverts$z == min(sverts$z), ]
            if (nrow(sverts_bot) == 1) {
                bot_r[i, j] <- sverts_bot$r[1]
                bot_g[i, j] <- sverts_bot$g[1]
                bot_b[i, j] <- sverts_bot$b[1]
            }
            
            # ok, now split up the points into those above or below
            # median value, which should split the verts in half
            # first check we have more than 1 point, this doesnt make sense with <2 row
            if (nrow (sverts) > 1) {
                med <- median (sverts$z)
                sverts_upper <- sverts[sverts$z > med, ]
                sverts_lower <- sverts[sverts$z <= med, ]
                
                # calc the metrics
                mean_r_upper[i, j] <- mean(sverts_upper$r)
                mean_g_upper[i, j] <- mean(sverts_upper$g)
                mean_b_upper[i, j] <- mean(sverts_upper$b)
                mean_r_lower[i, j] <- mean(sverts_lower$r)
                mean_g_lower[i, j] <- mean(sverts_lower$g)
                mean_b_lower[i, j] <- mean(sverts_lower$b)
            }
        }
    }
    print (paste ('row finished! :', i))
}
print ('all done!')

####################################################################
# FINAL CALCS AND PLOTS

# assemble rgbs
mean_rgb <- brick (c(mean_r, mean_g, mean_b))
top_rgb <- brick (c(top_r, top_g, top_b))
bot_rgb <- brick (c(bot_r, bot_g, bot_b))
mean_rgb_upper <- brick (c(mean_r_upper, mean_g_upper, mean_b_upper))
mean_rgb_lower <- brick (c(mean_r_lower, mean_g_lower, mean_b_lower))

# do range for proxy of local height
range_z <- max_z - min_z

# with these data you could go do some classification stuff that I haven't
# played around with here.

# make some jpegs (you could write geotiffs with writeRaster)
setwd ('..//jpegs//')
plot_to_file <- function (ras, filename) {
    # function to plot data to files, ras is raster (or brick if rgb)
    jpeg (filename)
    if (nlayers(ras) == 3) {
        plotRGB (ras)
    } else {
        plot (ras)
    }
    dev.off()
}

plot_to_file (mean_rgb, 'mean_rgb.jpg')
plot_to_file (top_rgb, 'top_rgb.jpg')
plot_to_file (bot_rgb, 'bot_rgb.jpg')
plot_to_file (mean_rgb_upper, 'mean_rgb_upper.jpg')
plot_to_file (mean_rgb_lower, 'mean_rgb_lower.jpg')
plot_to_file (mean_z, 'mean_z.jpg')
plot_to_file (max_z, 'max_z.jpg')
plot_to_file (min_z, 'min_z.jpg')
plot_to_file (range_z, 'range_z.jpg')
plot_to_file (numverts, 'numverts.jpg')


















