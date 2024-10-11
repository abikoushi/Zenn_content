#5*5
A <- matrix(0,24,24)
i <- 1
i+1
5*i + i

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

