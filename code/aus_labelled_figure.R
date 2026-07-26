# bring in all the bits for our map of Australian landmarks
# (code for the map itself is at the bottom but I'm calling this script from ecdfs.R)
library(RColorBrewer)
library(raster)
library(sf)
library(terra)

aus_simple <- st_read('data/aus_simplified.shp')

# grab an elevation raster :)
# library(elevatr)
# aus_sf = st_as_sf(aus_simple)
# aus_elevation = get_elev_raster(aus_sf, z=1) %>%
#   crop(aus_mask) %>%
#   projectRaster(aus_mask) %>%
#   raster::mask(aus_mask)

# murray darling basin boundary
# https://data.gov.au/data/dataset/murray-darling-basin-boundary/resource/8a6d889d-723b-492d-8c12-b8b0d1ba4b5a
mdb_shp <- shapefile("data/australian_landmarks/murray_darling_basin_boundary/mdb_boundary.shp")

# rivers?
# https://data.apps.fao.org/catalog/iso/6a53d768-1e20-46ea-92a8-c4040286057d
# predictably, this is huge
# simplify by strahler code? and write into a simplified shp so I never have to do this again?
rivers <- st_read("data/australian_landmarks/rivers_australia_37252/rivers_australia_37252.shp") %>%
  filter(grepl("Australia", MAJ_NAME) | MAJ_NAME == "Murray - Darling") %>%
  filter(Strahler > 3)


# clim_zones <- read.delim("data/australian_landmarks/kpn/kpnall.txt", sep=" ", skip=6) %>%
#   as.matrix() %>%
#   raster(xmn=112, ymn=-44, xmx=1681*0.025 + 112, ymx=1361*0.025 + -44) %>%
#   raster::mask(aus_simple)
# clim_zones[clim_zones == -9999] = NA
# # cut by major classes defined in the "all classes" README (which doesn't match the grouped READMEEEEE)
# values(clim_zones)[!is.na(values(clim_zones))] <- cut(values(clim_zones)[!is.na(values(clim_zones))],
#                                                       breaks=c(0,10,16,25,34,38,45))
# 
# bom_clim_pal = c("#63CAF3","#FFE493","#F79760","#D8E9B4","#8BAF3F","#B69535")
# clim_pal <- c("#DEC3D2","#F3E4B8","#F79760","#63CAF3","#14C061","#B69535")
# clim_pal <- c("#bebada","#ffffb3","#fdb462","#fb8072","#b3de69","#8dd3c7")

transf_seasrain <- function(from){
  c(1,1,1,2,3,3,3,2,4,4,4,2,5,5,5,2,6,6,6,2)[from]
}

clim_zones <- read.delim("data/australian_landmarks/seasb/seasrain.txt", sep=" ", skip=6) %>%
  as.matrix() %>%
  raster(xmn=112, ymn=-44, xmx=1681*0.025 + 112, ymx=1361*0.025 + -44)
clim_zones[clim_zones == -9999] = NA
# cut by major classes defined in the "all classes" README (which doesn't match the grouped READMEEEEE)
# values(clim_zones)[!is.na(values(clim_zones))] <- cut(values(clim_zones)[!is.na(values(clim_zones))],
#                                                       breaks=c(0,3,4,8,12,16,20))
values(clim_zones)[!is.na(values(clim_zones))] <- transf_seasrain(values(clim_zones)[!is.na(values(clim_zones))])

#plot(clim_zones, col=viridis(6)[sample(1:6)])
#plot(st_geometry(aus_simple), add=TRUE, border="white", lwd=2)

# attempt to de-blob:
w1 <- focalWeight(clim_zones, 1)
w1[w1!=0] = 1
w2 <- focalWeight(clim_zones, 0.3)
w2[w2!=0] = 1
clim_zones <- clim_zones %>%
  rast() %>%
  focal(w1, max, na.policy="only", na.rm=TRUE) %>% 
  # after much consternation, this is the only one I've found to fill in ugly NSW/SA hole
  disagg(fact = 2, method = "near") %>% # If I make this any higher, things get rather jaggedy
  focal(w1, median, na.rm=TRUE) %>% # why???
  disagg(fact = 5, method="near") %>%
  mask(aus_simple)

clim_pal <- c("#fb8072","#ffffb3","#fdb462","#b3de69","#8dd3c7", "#bebada")
#plot(clim_zones, col=clim_pal)


# city coords
city_coords <- read.csv("data/australian_landmarks/au.csv")
city_coords <- city_coords[city_coords$capital == "admin" | city_coords$capital == "primary",]
city_coords$labelx <- c(152, 148, 153.5, 115, 135, 151.5, 150, 127)
city_coords$labely <- c(-34, -39, -23, -32, -36, -36, -42, -12)
city_coords$pos <- c(4, 4, 3, 2, 2, 4, 4, 2)
city_coords$city[city_coords$city == "Canberra"] = "Canberra\n(ACT)"
city_coords <- rbind(city_coords, c("Tiwi Islands", -11.6347, 130.9228, 
                                    "Australia", "AU", "Northern Territory", 
                                    NA, NA, NA, 127,-11,2))
city_coords$labelx <- as.numeric(city_coords$labelx)
city_coords$labely <- as.numeric(city_coords$labely)
city_coords$lat <- as.numeric(city_coords$lat)
city_coords$lng <- as.numeric(city_coords$lng)
city_coords$pos <- as.numeric(city_coords$pos)

# state labels
state_coords <- data.frame("letters" = c("NSW", "NT", "QLD", "SA", "TAS", "VIC", "WA"),
                           "lon" = c(145,133.5,147,135,146.7,143,122),
                           "lat" = c(-32,-17,-24,-30,-42,-37,-25))

dark_purp <- colorRampPalette(c(brewer.pal(9,"Purples")[9], "black"))(10)[5]
purps = colorRampPalette(brewer.pal(9,"Purples"))(12)
blues = colorRampPalette(brewer.pal(9, "Blues")[2:9])(7)
#"#DEEBF7" "#BFD8EC" "#8DC0DD" "#56A0CE" "#2C7BBA" "#0C56A0" "#08306B"
blues = colorRampPalette(c("#dee1ff","#0317fc"))(5) # make the top blue more vibrant?
blues = colorRampPalette(c("#BDE2FF", "#007FDF"))(5)

curved_text = function(text, xlim, ylim, rotlim,
                       func, # get ycoords from here ..
                       dve, # get char rotations from here ..
                       inve, # get xcoords from here ..
                       func_over=c(0,1)){
  dummyx = seq(func_over[1], func_over[2], length.out=nchar(text))
  #return(dummyx)
  xcoords = inve(dummyx)
  xcoords = (xcoords - min(xcoords)) / (min(xcoords) - max(xcoords)) * -1
  ycoords = func(dummyx)
  ycoords = (ycoords - min(ycoords)) / (max(ycoords) - min(ycoords))
  srt = dve(dummyx)
  srt = (srt - min(srt)) / (max(srt) - min(srt))
  xcoords = xcoords * (xlim[2] - xlim[1]) + xlim[1]
  ycoords = ycoords * (ylim[2] - ylim[1]) + ylim[1]
  srt = srt * (rotlim[2] - rotlim[1]) + rotlim[1]
  
  return(data.frame(xcoords=xcoords,
                    ycoords=ycoords,
                    srt=srt,
                    char=strsplit(text, "")[[1]]))
}

# gdr = curved_text("Great Dividing Range", c(146,151), c(-37, -24), c(35, 135), 
#                   function(x){x}, function(x){x}, sin, func_over=c(1/8*pi, 6/8*pi))
# gdr = curved_text("Great Dividing Range", c(146,151), c(-37, -24), c(35, 135), 
#                   function(x){}, function(x){x}, function(x){sqrt(1 - x**2)}, func_over=c(-1,0.7))
# gdr = curved_text("Great Dividing Range", c(146,152), c(-37, -27), c(0, 90), 
#                   exp, function(x){log(x+0.1)}, function(x){log(x+0.1)})
gdr = curved_text("Great Dividing Range", c(146,151), c(-37, -24), c(0, 135), 
                  function(x){sin(x)},
                  function(x){x},
                  function(x){cos(x)},
                  func_over=c(-pi/2,pi/4))

# murray = curved_text("Murray River", c(140, 144), c(-36.2, -34.7), c(-45, 45),
#                      sin, cos, function(x){log(x+1)}, func_over=c(pi/4, 7/8*pi))
murray = curved_text("Murray River", c(140.2, 144.5), c(-36.4, -34.9), c(-45, 22),
                     sin, cos, function(x){log(x+1)}, func_over=c(pi*3/8, 7/8*pi))

# plot(aus_mask, xlim=c(140, 152), ylim=c(-40, -30))
# plot(aus_simple, add=TRUE)
# sapply(1:nrow(murray), function(x){
#   text(murray$xcoords[x], murray$ycoords[x], murray$char[x], srt=murray$srt[x])
# })


# fold everything into functions :)

plot_rivers <- function(lwd_fac=1){
  plot(st_geometry(rivers), col=blues[rivers$Strahler - min(rivers$Strahler) + 2], add=TRUE, 
       lwd=rivers$Strahler/min(rivers$Strahler)*lwd_fac)
}


plot_state_labels <- function(cex=1.3, alph=0.4){
  text(state_coords$lon, state_coords$lat, state_coords$letters, 
       col=alpha(dark_purp, alph), font=2, cex=cex)
}


plot_cities <- function(cex=1){
  sapply(1:nrow(city_coords), function(x){
    lines(c(city_coords$lng[x], city_coords$labelx[x]),
          c(city_coords$lat[x], city_coords$labely[x]))
    text(city_coords$labelx[x], city_coords$labely[x], city_coords$city[x],
         pos=city_coords$pos[x], cex=cex)
  })
  return(0)
}


plot_great_dividing_range <- function(cex=1){
  sapply(1:nrow(gdr), function(x){
    text(gdr$xcoords[x], gdr$ycoords[x], gdr$char[x], srt=gdr$srt[x], cex=cex)
  })
}


plot_murray <- function(){
  sapply(1:nrow(murray), function(x){
    text(murray$xcoords[x], murray$ycoords[x], murray$char[x], srt=murray$srt[x])
  })
}

#"#C23677"
plot_mdb <- function(col="#DE9DBB", alph=0.1, lwd=3, points=FALSE){
  if (points == TRUE){
    pts <- mdb_shp %>%
      st_as_sf() %>%
      st_transform(3112) %>%
      st_sample(size=500, type="hexagonal") %>%
      st_transform("+proj=longlat +datum=WGS84")
    plot(pts, add=TRUE, col=alpha(col, alph), cex=0.5, lwd=lwd)
    legend(x=127.5, y=-37.5, "Murray-Darling\nBasin", pch=1, 
           col=col, bty="n", pt.lwd=lwd+1, pt.cex=1.2)
  } else {
    plot(mdb_shp, add=TRUE, border=col, lwd=lwd, col=alpha(col, alph))
    legend(x=126.5, y=-38, "Murray-Darling Basin", pch=22, pt.bg=alpha(col, alph), 
           col=col, bty="n", pt.lwd=lwd+1, pt.cex=2.5)
  }
}


plot_trop_cap <- function(cex=1){
  lines(c(112, 155), rep(-23.4394, 2), lty=2)
  text(123, -23, "Tropic of Capricorn", font=3, cex=cex)
}

library("geosphere")
distHaversine(c(128,-25), c(133, -25))/1000

{png("~/Desktop/jev/pathogen_paper/figures/aus_labelled.png",
    height=1800, width=2400, pointsize=40)
  par(bty="n", mar=c(0.2,0.1,0.1,2.1))
  plot(raster(clim_zones), col=alpha(clim_pal, 0.4), legend=FALSE, xaxt="n", yaxt="n")
  # legend(117,-37, rev(c("Temperate", "Grassland", "Desert", "Subtropical", "Tropical", "Equatorial")),
  #        fill=rev(clim_pal), title="Koppen classification", bty="n")
  legend(117,-37, rev(c("Summer Dominant", "Arid", "Summer", "Uniform", "Winter", "Winter Dominant")),
         fill=rev(clim_pal), title="Long-term Rainfall Pattern", bty="n")
  lines(c(116.8, 136.5, 136.5, 116.8, 116.8), c(-37, -37, -44, -44, -37))
  axis(1, at=c(129, 131.5, 134), labels=c(0, 250, 500), line=-2.3) # distHaversine(c(128,-25), c(133, -25))/1000
  text(131.5, -41, "km")
# plot(aus_elevation, breaks=c(seq(0, maxValue(aus_elevation), 200), maxValue(aus_elevation)), col=purps,
#      legend=FALSE,
#      xaxt="n", yaxt="n")
# plot(aus_elevation, col=purps[1:7], legend.only=TRUE, border="grey80", lwd=2,
#      legend.args=list(text="Elevation", side=2, line=0.5), legend.width=1.2,
#      axis.args=list(at=c(minValue(aus_elevation) + 100, maxValue(aus_elevation) - 100), c("Low", "High")),
#      bty="o")
plot_mdb("#D780A9", points = TRUE, alph=0.7)
plot(st_geometry(aus_simple), add=TRUE, border="grey50", lwd=3)
plot_rivers(lwd_fac=1.5)
plot_state_labels()
plot_great_dividing_range()
plot_murray()
plot_cities()
text(143.5, -16, "Cape\nYork", adj=c(0.5, 0.5))
dev.off()}










