library(Matrix)
library(dplyr)
library(ggplot2)
library(gganimate)

make_adjmat_2d <-function(nr, nc){
  nn = as.integer(nr*nc)
  A <- Matrix::spMatrix(nrow = nn, ncol = nn)
  edges = nr*(0:nc) 
  for(i in 1L:nn){
    if(i-nr > 0){ #left bound
      A[i,i-nr] <- 1
    }
    if( i %% nr != 1L ){ #top bound
      A[i,i-1] <- 1
    }
    if((i+nr) <= nn){ #right bound
      A[i,i+nr] <- 1
    }
    if( (i %% nr != 0L) & i+1L <= nn ){ #bottom bound
      A[i,i+1] <- 1
    }
  }
  return(A)
}


difit <- function(coef, u, nt){
  U <- matrix(0, nt+1, length(u))
  U[1,] <- u
  for(i in 1:nt){
    u2 <- coef%*%u
    U[i+1,] <- as.vector(u2)
    u <- u2
  }
  return(U)
}

makecoef_exp <- function(r, A){
  diag(nrow = nrow(A)) + r*(A-diag(rowSums(A)))
}

makecoef_imp <- function(r, A){
  solve(diag(nrow = nrow(A)) - r*(A-diag(rowSums(A))))
}

A <- make_adjmat_2d(5,5)

jpeg("adjmat.jpg")
image(A)
dev.off()

r <- 0.02
#beta <- makecoef_exp(r, A)
beta <- makecoef_imp(r, A)

jpeg("coefmat.jpg")
image(beta)
dev.off()


uini <- numeric(length = nrow(A))
uini[1] <- 1
uini[8] <- 1

U <- difit(beta, uini, 24)

print(rowSums(U))


posid = expand.grid(y=1:5,x=1:5) %>% 
  mutate(pos=row_number())

dfU <- reshape2::melt(U) %>% 
  setNames(c("time","pos","density")) %>% 
  left_join(posid, by="pos")
dfU



ggplot(dfU,aes(x=x, y=density, group = time, colour=time))+
  geom_line()+
  facet_grid(rows = vars(y))+
  theme_bw(16)

ggsave("xoriented.png")

ggplot(dfU,aes(x=y, y=density, group = time, colour=time))+
  geom_line()+
  facet_grid(rows = vars(x))+
  theme_bw(16)

ggsave("yoriented.png")

ggplot(dfU,aes(x=x, y=y, fill=density))+
  geom_tile()+
  facet_wrap(~time)+
  scale_y_reverse()+
  scale_fill_gradient(low="white", high="black")+
  theme_bw(16)

  ggsave("tiles.png")


ggplot(dfU,aes(x=x, y=y, fill=density))+
  geom_tile()+
  scale_y_reverse()+
  scale_fill_gradient(low="white", high="black")+
  theme_bw(16)+
  transition_time(time)+labs(title =  'time: {frame_time}')

anim_save("tiles.gif")
