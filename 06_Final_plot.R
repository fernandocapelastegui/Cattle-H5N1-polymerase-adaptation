#Set parameters

font_size = 3
source("Code/colours.R")

#Import IQtree

tree <- read.tree(paste0("/Users/capelastegui.f/git/bovine_tree_figure/Data/04_trees/",tdate,"/",tdate,"_raw_concat_aligned_trimmed_clean.fasta.treefile"))

#####
#Make a generic tree we can interrogate with mutational profiles
#Load in the mutation data polymerase_geno so we have a tree object with data in

p<-ggtree(tree, branch.length = 'none', size=0.3) %<+% polymerase_geno  + 
  aes(color = concat)+
  geom_tippoint(aes(color = concat), size = 1) +  # Color only the tip labels
  scale_color_discrete("Mutation profile")+
  geom_text(aes(label=node), hjust=-.3, size = 2)
p

polymerase_heat_map$group <- as.factor(polymerase_heat_map$group)
polymerase_heat_map$branch_colour <- as.factor(polymerase_heat_map$branch_colour)

#Make a tree, with the highlighting based on the node numbers from tree t (heatmap) and the tree data enriched with mutation data (tree_data/polymerase_heat_map)
p<- ggtree(tree, branch.length = 'none', size=0.3) %<+% polymerase_heat_map+
    geom_tree(aes(color = branch_colour), show.legend = FALSE)+ 
  scale_color_manual(values = c("1" = "red", "2"="blue","0" = "black"))+
  geom_hilight(node=996, fill=colour_M631L, alpha=0.3)+
  geom_hilight(node=999, fill=colour_K497R, alpha=0.3)+
  geom_hilight(node=1645, fill=colour_I13V, alpha=0.9) +
  geom_hilight(node=192, fill=colour_D740N, alpha=1, extend=5) + #802
  geom_hilight(node=1251, fill=colour_D740N, alpha=0.9) +
  geom_hilight(node=292, fill=colour_D740N, alpha=1, extend=5) + #917
  geom_hilight(node=1039, fill=colour_Q591R, alpha=0.9)+
  geom_hilight(node=1653, fill=colour_E613K, alpha=0.9) +
  geom_hilight(node=1633, fill=colour_E677G, alpha=0.9) +
  
  geom_cladelab(node=996, label="PB2 M631L", align=TRUE, 
                geom='label', fill=colour_M631L,  barsize=0, hjust = 5, vjust = 20, fontsize = font_size, fontface="bold")+
  geom_cladelab(node=999, label="PA K497R", align=TRUE, 
                geom='label', fill=colour_K497R,  barsize=0, hjust =5.7, vjust = 18, fontsize = font_size, fontface="bold")+
  geom_cladelab(node=1645, label="PA I13V", align=TRUE, 
                geom='label', fill=colour_I13V, textcolour = "white", fontsize = font_size, vjust =-0.4,fontface="bold")+
  geom_cladelab(node=c(192), label="PB2 D740N", align=TRUE, 
                geom='label', fill=colour_D740N, fontsize = font_size, textcolour = "white", fontface="bold")+
  geom_cladelab(node=c(1251), label="PB2 D740N", align=TRUE, 
                geom='label', fill=colour_D740N, fontsize = font_size, textcolour = "white", fontface="bold", vjust = 0.2)+
  geom_cladelab(node=c(292), label="PB2 D740N", align=TRUE, 
                geom='label', fill=colour_D740N, fontsize = font_size, textcolour = "white", fontface="bold")+
  geom_cladelab(node=1039, label="PB2 Q591R", align=TRUE, 
                geom='label', fill=colour_Q591R, fontsize = font_size, fontface="bold")+
  geom_cladelab(node=1653, label="PA E613K", align=TRUE, 
                geom='label', fill=colour_E613K, fontsize = font_size, fontface="bold", vjust =0.2)+
  geom_cladelab(node=1633, label="PB2 E677G", align=TRUE, 
                geom='label', fill=colour_E677G, fontsize = font_size, fontface="bold")+
  geom_tiplab(aes(label = ifelse(label == c(reference_name), reference_name_short, "")),
              color = "red", fontface = "bold", size = 2.5, vjust =0.5, hjust = -0.71) +
  geom_tippoint(aes(subset = (label == c(reference_name))), color = "red", size = 3)+
  geom_tiplab(aes(label = ifelse(label == c(texas37_name), texas37_name_short, "")),
              color = "blue", fontface = "bold", size = 2.5, hjust =-0.965, vjust=-0.5) +
  geom_tippoint(aes(subset = (label == c(texas37_name))), color = "blue", size = 2.2, shape = 17)+
  xlim(-50,400)+
  ylim(-50,930)+
  geom_strip("EPI_ISL_18698985_A/gyrfalcon/Idaho/23_034082_001/2023_2023_10_27", 
             "EPI_ISL_19014396_A/Wild_Bird/Wyoming/24_003692_001/2024_2024_01_25", 
  barsize=2, color='orange', 
             label = "US avian sequences\nSeptember 2023 - March 2024", 
  offset = 22,
  offset.text=2)

 p
#Save out the tree data to idenitify clades:
tree_data<-p[["data"]]

#Add in picture symbols for goose and cow
phylopic_info <- data.frame(node = c(774,775,998),
                            phylopic = c("036b96de-4bca-408e-adb6-2154fcd724ef",
                                         "b677ec7b-a2ef-46be-9d78-a997a88e1a9c",
                                         "8a2d5863-cd6c-4158-a164-d62e789c2210"),
                            species = c("human", "goose", "cow"),
                            colour = c("blue", "red", "black"))

phylopic_info_goose <- phylopic_info %>% filter(species=="goose") %>% 
  mutate(tip.label = reference_name)

p1<-p %<+% phylopic_info_goose + geom_tiplab(aes(image=phylopic), 
                                             geom="phylopic", 
                                             alpha=1, 
                                             color = "red", 
                                             nudge_x = 105,
                                             #nudge_y = -10,
                                             size=0.05)
p1


phylopic_info_human<- phylopic_info %>% filter(species=="human") %>% 
  rename(species3 = species,
         phylopic3 = phylopic)

p2 <- p1 %<+% phylopic_info_human + geom_tiplab(aes(image=phylopic3), 
                                                geom="phylopic", 
                                                alpha=1, 
                                                color = "blue", 
                                                nudge_x = 95,
                                                #nudge_y = -10,
                                                size=0.025)
p2


phylopic_info_cow<- phylopic_info %>% filter(species=="cow") %>% 
  rename(species2 = species,
         phylopic2 = phylopic)



p3<- p2 + geom_cladelab(data = phylopic_info_cow,
                        mapping = aes(image = phylopic2, node = node, label = species2),
                        imagecolor="black",
                        geom = "phylopic", 
                        offset.text = 10, 
                        show.legend = FALSE, 
                        offset = 22,
                        imagesize = 0.075,
                        barsize = 2
)
p3

  #p3<- p3+theme(plot.margin = margin(40, 40, 40, 40))

p3
#Save out
#ggsave(paste0("/Users/capelastegui.f/git/bovine_tree_figure/Outputs/",tdate,"/",tdate,"bigger_plot.png"),dpi = 600)
ggsave(paste0("/Users/capelastegui.f/git/bovine_tree_figure/Outputs/",tdate,"/",tdate,"_Final_plot.pdf"),
        width = 50, height = 30, units = "cm", limitsize = FALSE)
