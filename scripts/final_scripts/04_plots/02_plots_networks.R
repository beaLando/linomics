# LIBRARIES
library(ggplot2)
library(tidyverse)
library(ape)
library(ggtree)
library(tanggle)

# Read your extended Newick (from PhyloNet output) - starting tree provided
networks.ML.none <- data.frame(h = c(0, 1, 2, 3, 4, 5, 6),
                               loglik = c(-532.5543974361876,
                                          -519.5295512001662,
                                          -520.4758929101513,
                                          -513.6384766702216,
                                          NA,
                                          NA,
                                          NA),
                               hybrid_node_set = "None")

networks.ML.lb <- data.frame(h = c(0, 1, 2, 3),
                           loglik = c(-532.5543974361876,
                                      -524.566708409653,
                                      -522.590138487594,
                                      -521.7475171480611
                                      ),
                           hybrid_node_set = "L. bienne") 

networks.ML.lblh <- data.frame(h = c(0, 1, 2, 3),
                             loglik = c(-532.5543974361876,
                                        -522.4475838926337,
                                        -522.658903539097,
                                        -522.0950911365452
                             ),
                             hybrid_node_set = "L. bienne and L. hologynum") 

networks.ML <- bind_rows(networks.ML.none,
                         networks.ML.lb) %>%
  bind_rows(.,
            networks.ML.lblh) %>%
  as.data.frame()


networks.ML.lbonly.none <- data.frame(h = c(0, 1, 2, 3),
                               loglik = c(-319.2550949719137,
                                          -317.55568832919204,
                                          -314.57860220799296,
                                          -314.579466020863),
                               hybrid_node_set = "None")

networks.ML.lbonly.lb <- data.frame(h = c(0, 1, 2, 3),
                             loglik = c(-319.2550949719137,
                                        -316.08160817193084,
                                        -313.50232505192827,
                                        -313.5408410924158
                             ),
                             hybrid_node_set = "L. bienne") 

networks.ML.lbonly <- bind_rows(networks.ML.lbonly.none,
                         networks.ML.lbonly.lb) %>%
  as.data.frame()


# PLOTS
## CRITERIA
ggplot(networks.ML,
       aes(x = as.factor(h),
           col = hybrid_node_set)) +
  geom_point(aes(y = loglik)) +
  geom_line(aes(y = loglik, group = hybrid_node_set)) +
  scale_colour_manual(values = c("black", "red", "blue")) +
  labs(x = "Maximum Number of Hybrid Nodes (h)", y = "Log Likelihood", col = "Hybrid Node Set To") +
  theme_classic(base_size = 18)

ggplot(networks.ML.lbonly,
       aes(x = as.factor(h),
           col = hybrid_node_set)) +
  geom_point(aes(y = loglik)) +
  geom_line(aes(y = loglik, group = hybrid_node_set)) +
  scale_colour_manual(values = c("red", "black")) +
  labs(x = "Maximum Number of Hybrid Nodes (h)", y = "Log Likelihood", col = "Hybrid Node Set To") +
  theme_classic(base_size = 18)


# PLOT
## NETWORKs - W/ LB only
### Networks ML - LB as hyb
trML1_lb <- read.evonet(text="(hirsutum:0.12690773023966948,((bienne:1.0902006095730272)#H1:1.8724269407494951::0.050808538302931416,(narbonense:0.022391463211522272,(#H1:2.3734711124382297::0.9491914616970686,(grandiflorum:2.126756343578251,decumbens:7.371305173148464):2.0732583843552446):0.15282217871639758):4.913000167045098):2.626155661917847);")

ggevonet(trML1_lb, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML1b_lb <- read.evonet(text="(hirsutum:0.24631933296679265,(narbonense:0.32101000233604277,((bienne:2.995739592829998)#H1:0.21176525788609624::0.7469064408194991,(#H1:0.20641110044257371::0.2530935591805009,(grandiflorum:0.010090039570159914,decumbens:3.747682004099284E-4):2.3181690854699517):0.04568291440270951):0.101163336399686):8.29617966132559);")
#-317.5273378506653
ggevonet(trML1b_lb, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML2_lb <- read.evonet(text="(hirsutum:0.09550543415927346,(((bienne:0.17768879452145242)#H2:0.916633414416105::0.03146842400354932)#H1:1.5250760996086383::0.6754150706018403,(narbonense:0.043721240616240775,((decumbens:2.5296270732226764,grandiflorum:0.3807634014267618):2.321220293756171,(#H2:4.530240370173104::0.9685315759964507,#H1:1.934670706638734::0.32458492939815975):1.196582682223942):0.15735632523106577):1.848114439340342):6.2638035885686625);")

ggevonet(trML2_lb, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML2b_lb <- read.evonet(text="(hirsutum:0.13778764445329944,((bienne:1.143358990489759)#H1:0.6776406575113579::0.034526188711206274,(narbonense:1.8711374828142109,((#H1:2.8961027672530775::0.9654738112887937)#H2:2.2246648433438283::0.5315546325682967,(#H2:1.2940081871963063::0.46844536743170323,(decumbens:0.6382740611092623,grandiflorum:0.15734896904755435):2.6945602445222194):0.09459574058826987):0.10473207446787859):13.624639655943074):9.4623259539236);")

ggevonet(trML2b_lb, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

### Networks ML - All (no hyb set)
trML0_all <- read.tree(text="(hirsutum:0.1579930565744506,(narbonense:1.0464525097746953,(bienne:2.7133935068022956,(decumbens:3.675917020921333,grandiflorum:7.228650119436454):2.5497172451091403):0.20141342693150138):7.644079301642109);")

ggtree(trML0_all) + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML1_all <- read.evonet(text="(hirsutum:2.588124154455026,((narbonense:4.352561007397007)#H1:4.160980204388217::0.3593330663117412,(((grandiflorum:1.0626549158466772,decumbens:1.791714742640742):2.428424933176867,bienne:2.7742104604406226):0.03726829196603053,#H1:0.28783354515068804::0.6406669336882588):0.4820134118361771):6.815312344154125);")

ggevonet(trML1_all, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML2_all <- read.evonet(text="(hirsutum:2.4300155994641384,((bienne:1.1931621530863765)#H1:0.32600995962292667::0.06326586798729772,(narbonense:2.168578258260648,((grandiflorum:2.836777890759998,decumbens:4.98235705267143):2.384356714118305,#H1:3.697433282569384::0.9367341320127023):0.10940978335859466):0.7171571465340293):6.434496922363379);")

ggevonet(trML2_all, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML2b_all <- read.evonet(text="(hirsutum:1.372319448718141,(((decumbens:1.1406570948910588,grandiflorum:0.46205846247246574):0.90165116060997)#H1:1.845205762323045::0.11437062717380829,(((bienne:2.9724654590006874,#H1:5.9256672745653125::0.8856293728261917):3.2626138617583567E-4)#H2:1.0::0.5,(narbonense:1.4250979554699779,#H2:0.144002305452711::0.5):0.23141051740462837):0.3212076267057903):7.002601395276301);")
#LL=-316.47054760049207
ggevonet(trML2b_all, layout = "rectangular") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()


## NETWORKs - W/ also LH
### Networks ML - LH and LB as hyb
trML1_lblh <- read.evonet(text="(hirsutum:4.9982550456242905,(((bienne:2.03538702365111,hologynum:2.0):0.11861277653925006)#H1:0.0011282902474240575::0.13724077376502597,(narbonense:5.873857292690451,((grandiflorum:4.581726544523024,decumbens:1.8638321197544965):2.310281469314684,#H1:3.0874782233697355::0.862759226234974):0.16567969268012978):0.32286659259357897):7.334855888384153);")

ggevonet(trML1_lblh, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML2_lblh <- read.evonet(text="(hirsutum:1.6253825207508958,((hologynum:0.012794140139675975)#H1:3.771199613526808::0.12722246635284862,((bienne:1.887091233472067)#H2:0.028371092805784696::0.1271893976255305,(narbonense:0.3456954393864071,((grandiflorum:11.73087110344752,decumbens:1.1297762859972242):2.045531138433374,(#H1:0.04443439370085417::0.8727775336471514,#H2:0.9614097669217547::0.8728106023744695):4.25679752292667):0.16368635225842584):0.09135583748726252):0.0968076827031655):11.987150077073215);")

ggevonet(trML2_lblh, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

### Networks ML - LB as hyb
trML1_lb <- read.evonet(text="(hirsutum:0.07229508902521778,((bienne:1.6005756698122018)#H1:1.380363203681051::0.04693742562450498,(narbonense:1.1596986520298016,((grandiflorum:2.9747083854094507,decumbens:0.09668540595153247):2.5016913449064933,(hologynum:1.691994079191224,#H1:0.8941659769943799::0.953062574375495):1.3247447856852417):0.14034449967394652):0.45687283308204):9.725075289251683);")

ggevonet(trML1_lb, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML2_lb <- read.evonet(text="(hirsutum:2.298459886738372,(((bienne:0.4488777929181677)#H2:0.3580411752573226::0.023392225954576684)#H1:0.9018134123110118::0.5,(narbonense:4.908916791087957,(#H1:1.0::0.5,((grandiflorum:0.4087113353422961,decumbens:1.042524786361489):2.321503760621748,(hologynum:4.939719271587901,#H2:2.4691616358185935::0.9766077740454233):1.2615769775595103):0.06460505793921295):0.04519819897983624):5.719441676328853):2.4199481361551887);")

ggevonet(trML2_lb, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

### Networks ML - All (no hyb set)
trML0_all <- read.tree(text="(hirsutum:1.2708194908838562,((hologynum:1.8791583480219205,bienne:2.1994478993914566):1.2314849623004707,((grandiflorum:7.305254376944676,decumbens:0.8840095857049515):2.495623048837901,narbonense:0.015012221529676423):0.016661232410812604):14.013743971084018);")

ggtree(trML0_all) + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML1_all <- read.evonet(text="(hirsutum:0.20447079749953734,(((narbonense:1.287912354754834,(grandiflorum:2.5235508830008917,decumbens:2.3238494259098315):2.432656069420193):0.13244243353943294)#H1:0.1275340765411558::0.4821231359770699,((bienne:2.046348565092514,hologynum:2.367385795225082):0.8623273316640256,#H1:0.056533346666477094::0.5178768640229301):2.33574355013411):7.675068743922907);")

ggevonet(trML1_all, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()

trML2_all <- read.evonet(text="(hirsutum:1.1107199810460964,(((narbonense:1.6382161348588427,(grandiflorum:2.874365681644031,decumbens:1.3114409021548528):2.6354349896362645):0.03594624034465137)#H1:0.09691932667190419::0.6861638269460952,(#H1:1.8200636941668642::0.3138361730539047,((bienne:2.2748397446812967)#H2:0.9056074644971257::0.5203241693260997,(#H2:0.924835328776565::0.4796758306739003,hologynum:0.0769353665472425):0.15635998806955326):0.5083707145267352):1.7507527297098708):8.322399961690717);")

ggevonet(trML2_all, layout = "slanted") + 
  geom_nodepoint() +
  geom_tiplab() + 
  geom_nodelab()
