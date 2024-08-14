library(readr)
library(sf)
library(tidyverse)
library(tmap)
library(lwgeom)
library(zipangu)
# library(gganimate)
library(animation)
# library(ggraph)
# require(tidygraph)
library(arcdiagram)
#devtools::install_github('gastonstat/arcdiagram')

un_graphe <- rbind(
  c("fromage", "pain"),
  c("pain", "vin"),
  c("vin", "biere"),
  c("cidre", "biere"),
  c("foie", "fromage"),
  c("pain", "foie"))

arcplot(un_graphe)
arcplot(un_graphe, horizontal=FALSE)
arcplot(un_graphe, above = c(1, 3, 5))
arcplot(un_graphe, sorted=TRUE)

path <- "/Users/koabe/data/gm-jpn-all_u_2_2/"
gm_jpn <- read_sf(paste0(path, "polbnda_jpn.shp"))

class(gm_jpn)

gm_jpn <- group_by(gm_jpn,nam) %>%
  summarize()
head(gm_jpn)
#plot(gm_jpn)

ggplot(data=gm_jpn) +
  geom_sf() + 
  coord_sf(crs=sf::st_crs("EPSG:3857"),
           default_crs = sf::st_crs("EPSG:4326")) + 
  theme_bw()


jpname <- read_csv(paste0(path,"data-map-JIS.txt"), col_names = c("num","name"))
jpmap <- read_csv(paste0(path,"data-map-neighbor.txt"))
jpname
range(jpmap$From)
range(jpmap$To)

A <- matrix(0,46,46)
A[cbind(jpmap$From,jpmap$To)] <- 1
A[cbind(jpmap$To, jpmap$From)] <- 1
image(A)
L <- diag(rowSums(A)) - A
image(L)
L
nt <- 19
set.seed(1234)
u <- rbinom(46,1,0.2)

difit <- function(L,u,nt,r){
  U <- matrix(0, nt+1, length(u))
  U[1,] <- u
  r <- 0.05
  for(i in 1:nt){
    u2 <- u - r*L%*%u
    U[i+1,] <- u2
    u <- u2
  }
  return(U)
}

dfU <- reshape2::melt(U) %>% 
  setNames(c("time","num","density")) %>% 
  mutate(jis_code = if_else(num<10,paste0("0",num),as.character(num))) %>% 
  left_join(jpnprefs, by="jis_code") %>% 
  mutate(nam=gsub("-k"," K", prefecture)) %>% 
  left_join(gm_jpn, by="nam")

dfU
ggplot(dfU,aes(x=reorder(prefecture_kanji, num), y=density, group=time, colour=time))+
  geom_line()+
  scale_colour_viridis_c()+
  theme(axis.text.x = element_text(family="Osaka", angle=90, hjust=1))+
  labs(x="")
ggsave("diff_jp.png")

head(dfU)


saveGIF({
  for(i in 1:20){
    p1 <- ggplot(dfU[dfU$time==i,], aes(geometry=geometry, fill=density)) +
      geom_sf() + 
      theme_minimal() +
      theme(axis.text = element_blank(), axis.ticks = element_blank())+
      facet_wrap(~time, labeller = label_both)+
      scale_fill_gradient(low="white", high="royalblue", limits = c(0,1))
    print(p1)
  }
})

