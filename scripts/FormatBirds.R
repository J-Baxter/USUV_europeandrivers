# 
isBird <- function(x){
  bird_orders <- c(
    "passeriformes", "charadriiformes", "suliformes", "pelecaniformes","accipitriformes", 
    "falconiformes", "galliformes", "piciformes", "anseriformes", "coraciiformes", "strigiformes", 
    "columbiformes", "gruiformes","procellariiformes", "ciconiiformes", "podicipediformes", 
    "apodiformes", "caprimulgiformes", "psittaciformes", "sphenisciformes", "casuariiformes",
    "phoenicopteriformes", "cuculiformes", "coraciformes", "trogoniformes", "gaviiformes", 
    "pteroclidiformes", "anseranseriformes", "phaethontiformes","opisthocomiformes")
  
  if(any(x %in% bird_orders)){
    x <- TRUE
  }else{
    x <- FALSE
  }
  
  return(x)
}


FormatAnseriformes <- function(x){
  anseriformes <- c(
    "whooper swan" = "cygnus cygnus|whooper swan",
    "mute swan" = "c[yi]gnus olor|m(ute){0,1}[ -]swan",
    "domestic goose sp. (domestic type)" = "(domestic|pomeranian|embden|sebastopol|american buff|african|rural) goose|anser anser domesticus",
    "mallard (domestic type)" = "(domestic|pekin|mule|runner|cascade|rural|mulard) duck",
    "barnacle goose" = "branta leucopsis|barnacle goose|^barnacle$",
    "mandarin duck" = "aix galericulata|man{0,1}darin duck",
    "duck sp." = "anatidae|^duck$|(wild|migratory) duck|^dk$|duck sp\\.",
    "greater/lesser white-fronted goose" = "^white[- ]fronted goose$",
    "goose sp." = "anatidae \\(goose sp\\.\\)|^goose$|wild (goose|geese)|^geese$|^go$",
    "greylag goose" = "anser anser|^(gr[ea]y {0,1}lag) goose$|grey go*",
    "green-winged teal" = "anas crecca|greenwing duck|(common|eurasian)[ -]teal|green[- ]{0,1}winged[- ]{0,1}teal$",
    "muscovy duck" = "cairina moschata|muscovy duck", #
    "mallard" = "anas platyrhynchos|mallard duck|^mal$|plath{0,1}yrh{0,1}ynchos|^mallard$",
    "swan sp." = "^cygnus$|^swan$",
    "greater white-fronted goose" = "anser albifrons|greater[- ]white[- ]fronted[- ]goose",
    "northern pintail" = "anas acuta|^pintail$|northern pintail",
    "canada goose" = "branta canadensis|canad[ae] goose",
    "whistling-duck sp." = "dendrocygna|whistling[- ]duck",
    "falcated duck" = "mareca falcata|^falcated (teal|duck)$",
    "brant" = "branta bernicla|^brant$",
    "egyptian goose" = "alopochen aegyptiaca|egyptian goose",
    "gadwall" = "mareca strepera|^gadwall*",
    "pink-footed goose" = "anser brachyrhynchus|pink[- ]footed goose",
    "black swan" = "cygnus atratus|(black|bk)[- ]swan",
    "tundra swan" = "cygnus columbianus|^(tundra|bewicks{0,1}) swan$",
    "ferruginous duck" = "aythya nyroca|ferruginous duck",
    "garganey" = "spatula querquedula|garganey",
    "taiga/tundra bean-goose" = "anser fabalis|anser serrirostris|bean[ -]goose|taiga/tundra bean-goose",
    "teal sp." = "^teal$",
    "eastern spot-billed duck" = "anas zonorhyncha|spot-billed duck|eastern spot-billed duck",
    "tufted duck" = "aythya fuligula|tufted duck|t dk",
    "eurasian/american wigeon" = "mareca(penelope|americana)|^wigeon$|(eur|european|eurasian) wig(eon){0,1}",
    "hawaiian goose" = "branta sandvicensis|(hawaiian|nene) goose",
    "waterfowl sp." = "anseriformes|waterfowl|wild waterbird",
    "swan goose" = "anser cygnoides|swan goose",
    "american wigeon" = "mareca americana|american wigeon",
    "northern shoveler" = "spatula clypeata|(northern ){0,1}shoveler",
    "bufflehead" = "bucephala albeola|bufflehead",
    "wood duck" = "aix sponsa|(wood|carolina) duck",
    "lesser scaup" = "aythya affinis|lesser scaup",
    "snow goose" = "anser caerulescens|snow goose",
    "blue-winged teal" = "spatula discors|blue[ -]winged[ -]teal",
    "red-breasted goose" = "branta ruficollis|red[- ]breasted goose",
    "ross's goose" = "anser rossii|ross'{0,1}s goose",
    "merganser sp." = "mergus|merganser",
    "dabbling duck sp." = "shelduck",
    "steamer-duck sp." = "steamer[- ] duck",
    "bar-headed goose" = "anser indicus|bar[- ]headed goose",
    "cape teal" = "anas capensis|cape teal",
    "coscoroba swan" = "coscoroba coscoroba|coscoroba swan",
    "trumpeter swan" = "cygnus buccinator|trumpeter swan",
    "cackling goose" = "branta hutchinsii|cackling goose",
    'brant' = '^br[ea]nt( goose){0,1}$|branta bernicla',
    'aythya sp.' = 'aythya sp\\.|por{0,1}chard$',
    'indian spot-billed duck' = 'indian spot[- ]billed duck|anas poecilor{0,1}hyncha|spotbill duck',
    'common goldeneye' = 'common golden {0,1}eye|bucephala clangula',
    'common eider' = "common eider|somateria mollissima|(st\\. cuthbert's|cuddy's) duck",
    'redhead' = 'redhead( duck){0,1}|aythya americana',
    'surf scoter' = 'surf sco{1,2}ter|melanitta perspicillata',
    'common scoter'='melanitta nigra|common scoter',
    'steamer-duck sp.'= "tachyeres( sp\\.){0,1}|steamer[- ]duck",
    'common shelduck' =	'tadorna tadorna|common shelduck'
    )

  for (i in 1:length(anseriformes)){
    if(any(grepl(anseriformes[[i]], x))){
      x <- names(anseriformes)[[i]]
    }
  }
  
  return(x)
}


FormatCharadriiformes <- function(x){
  charadriiformes <- c(
    "eurasian curlew" = "numenius arquata|eurasian curlew",
    "black-headed gull" = "chroicocephalus ridibundus|blac[kh][- ]headea{0,1}d gull|bl h gull",
    "eurasian oystercatcher" = "haematopus ostralegus|eurasian oystercatcher",
    "lapwing sp." = "vanellus sp\\.|lapwing sp.|lapwing",
    "red knot" = "calidris canutus|red knot|knot wader",
    "herring gull" = "larus argentatus|(herring|herrin) gull",
    "sanderling" = "calidris alba|sanderling",
    "oystercatcher sp." = "haematopus sp\\.|oystercatcher sp.|^oyster {0,1}catcher$",
    "gull sp." = "larus sp\\.|gull sp\\.|^c gull$|^sea {0,1}gull$|^gull\\d{0,2}$|^black[ -]backed[ -]gull",
    "brown-headed gull" = "chroicocephalus brunnicephalus|brown[- ]headed gull",
    "great skua" = "stercorarius skua|great skua",
    "curlew sp." = "numenius sp\\.|curlew sp\\.|^curlew$",
    "yellow-legged gull" = "larus michahellis|yellow-legged gull",
    "whiskered tern" = "chlidonias hybrida|whiskered tern",
    "ruddy turnstone" = "arenaria interpres|ruddy turnstone|^turnstone$",
    "gull/tern sp." = "laridae|ternidae|gull/tern sp\\.|^tern$",
    "caspian gull" = "larus cachinnans|caspian gull",
    "common gull" = "larus canus|common gull",
    "hartlaub's gull" = "chroicocephalus hartlaubii|hartlaub'{0,1}s gull",
    "kelp gull" = "larus dominicanus|kelp gull",
    "long-tailed jaeger" = "stercorarius longicaudus|long-tailed (jaeger|skua)",
    "little gull" = "hydrocoloeus minutus|little gull",
    "sandwich tern" = "thalasseus sandvicensis|sandwich tern",
    "great crested tern" = "thalasseus bergii|great crested tern|swift tern",
    "african oystercatcher" = "haematopus moquini|african (black ){0,1}oystercatcher",
    "arctic tern" = "sterna paradisaea|arc{0,1}tic tern",
    "ring-billed gull" = "larus delawarensis|ring-billed gull",
    "franklin's gull" = "larus pipixcan|franklin'{0,1}s gull",
    "common tern" = "sterna hirundo|common tern",
    "slaty-backed gull" = "larus schistisagus|slaty-backed gull",
    "black skimmer" = "rynchops niger|black skimmer",
    "glaucous-winged gull" = "larus glaucescens|glaucous-winged gull",
    "black-legged/red-legged kittiwake" = "rissa tridactyla|black-legged/red-legged kittiwake|^kittiwake$",
    "mediterranean gull" = "ichth{0,1}yaetus melanocephalus|mediterranean gull",
    "slender-billed gull" = "larus genei|slender-billed gull",
    "glaucous gull" = "larus hyperboreus|glaucous gull",
    "black-tailed gull"	= "larus crassirostris|black-tailed gull",
    "western gull" = "larus occidentalis|western gull",
    "far eastern curlew" = "numenius madagascariensis|*eastern curlew",
    "common murre" = "uria aalge|(common|foolish) guillemot|common murre",
    'great black-backed gull' = 'gr bk bd gull|(great[ -])black[ -]backed[ -]gull',
    "belcher's gull" = "larus belcheri|belcher'{0,1}s{0,1} gull",
    'lesser black-backed gull' = "lesser black[- ]backed gull|l[ -]bl[ -]ba[ -]gull",
    'razorbill' = "razorbill|lesser auk|alca torda",
    'large alcid sp.' = '^guillemot$',
    'black-bellied plover' = 	"pluvialis squatarola|(black[- ]belli(ed){0,1}|gr[ea]y) plover",
    'scolopacidae sp.' = 'sandpiper',
    'shorebird sp.' = 	'^shorebird$|charadriiformes sp\\.|^seabird$'
    
 
  )
  
  for (i in 1:length(charadriiformes)){
    if(any(grepl(charadriiformes[[i]], x))){
      x <- names(charadriiformes)[[i]]
    }
  }
  
  return(x)
}


FormatAccipitriformes <- function(x){
  accipitriformes <- c(
    "hawk sp." = "accipitridae sp.|hawk sp.|^hawk$",
    "common buzzard" = "buteo buteo|common buzzard|^buzzard$",
    "eastern buzzard" = "buteo japonicus|eastern buzzard",
    "eurasian goshawk" = "accipiter gentilis|(eurasian|northern) goshawk|^goshawk$",
    "white-tailed eagle" = "haliaeetus albicilla|(white-tailed|withe-tiled{0,1}) eagle",
    "golden eagle" = "aquila chrysaetos|golden eagle",
    "western marsh harrier" = "circus aeruginosus|western marsh harrier",
    "eagle sp." = "accipitridae \\(eagle sp\\.\\)|eagle sp.|sea eagle|^eagle$",
    "eurasian sparrowhawk" = "accipiter nisus|eurasian sparrowhawk|^sparrowhawk$",
    "cooper's hawk" = "accipiter cooperii|coopers'{0,1} {0,1}s{0,1} hawk",
    "bald eagle" = "haliaeetus leucocephalus|bald eagle",
    "red-tailed hawk" = "buteo jamaicensis|(red-tailed|rt)[- ]hawk",
    "african fish-eagle" = "haliaeetus vocifer|african fish[- ]eagle|^fish[- ]eagle$",
    #"old world vulture sp." = "accipitridae \\(old world vulture sp\\.\\)|old world vulture sp.",
    "variable hawk" = "buteo polyosoma|(variable|red[- ]backed)[- ]hawk",
    "bearded vulture" = "gypaetus barbatus|bearded vulture",
    "red-shouldered hawk" = "buteo lineatus|red-shouldered hawk",
    "osprey" = "pandion haliaetus|osprey",
    "harris's hawk" = "parabuteo unicinctus|harris'{0,1}s{0,1} hawk", 
    "rough-legged hawk" = "buteo lagopus|rough[- ]legged (hawk|buzzard)",
    "swainson's hawk" = "buteo swainsoni|swainson'{0,1}s{0,1} hawk"
  )
  
  for (i in 1:length(accipitriformes)){
    if(any(grepl(accipitriformes[[i]], x))){
      x <- names(accipitriformes)[[i]]
    }
  }
  
  return(x)
}


FormatPelecaniformes <- function(x){
  pelecaniformes <- c(
    "great egret" = "ardea alba|great egret",
    "grey heron" = "ardea cinerea|gr[ea]y heron",
    "white egret sp." = "egretta sp.|white egret sp.|^egret$",
    "great white pelican" = "pelecanus onocrotalus|great[- ]white[- ]pelican",
    "dalmatian pelican" = "pelecanus crispus|dalmatian pelican",
    "pelican sp." = "pelecanidae \\(pelican sp\\.\\)|pelican sp\\.|pelecanus|^pelican$",
    "eurasian spoonbill" = "platalea leucorodia|eurasian spoonbill",
    "ibis sp." = "threskiornithidae \\(ibis sp\\.\\)|ibis sp\\.|^ibis$",
    "african/malagasy sacred ibis" = "sacred ibis|threskiornis (aethiopicus|bernieri)",
    "heron sp." = "ardeidae \\(heron sp\\.\\)|heron sp\\.|^heron$",
    "black-headed heron" = "ardea melanocephala|black[- ]headed heron",
    "american white pelican" = "pelecanus erythrorhynchos|american (white ){0,1}pelican",
    "little egret" = "egretta garzetta|little egret",
    "great blue heron" = "ardea herodias|great blue heron",
    "brown pelican" = "pelecanus occidentalis|brown pelican",
    "snowy egret" = "egretta thula|snowy egret",
    'black-crowned night heron'	= 'nycticorax nycticorax|black[- ]crown(ed){0,1} night heron'
  )
  
  for (i in 1:length(pelecaniformes)){
    if(any(grepl(pelecaniformes[[i]], x))){
      x <- names(pelecaniformes)[[i]]
    }
  }
  
  return(x)
}


FormatGalliformes <- function(x){
  galliformes <-  c(
    "red junglefowl (domestic type)" = "gallus gallus|red junglefow|c[hk]i[ck][hk]{0,1}en|(lay(er|ing)|breeding) hen|broiler|poultry|^hen$|^rooster$|^layer$|^ch$|backyard bird",
    "wild turkey (domestic type)" = "meleagris gallopavo|wild turkey|turkey|^pavo$|domestic turkey",
    "helmeted guineafowl (domestic type)" = "numida meleagris|helmeted guineafowl \\(domestic type\\)",
    "pheasant sp." = "phasianidae|pheasant sp\\.|^pheasant$|^partridge$",
    "indian peafowl (domestic type)" = "pavo cristatus|indian peafowl|pea(fowl|cock)",
    #"old world quail sp." = "phasianidae (quail sp\\.\\)|old world quail sp.",
    #"new world quail sp." = "odontophoridae (quail sp\\.\\)|new world quail sp.",
    "lady amherst's pheasant" = "chrysolophus amherstiae|lady amherst'{0,1}s pheasant",
    "ring-necked pheasant" = "phasianus colchicus|(ring[- ]{0,1}neck(ed){0,1}|common) pheasant", 
    "crested guineafowl sp." = "guttera pucherani|crested guineafowl sp\\.",
    "red-legged partridge" = "alectoris rufa|red-legged partridge")

  for (i in 1:length(galliformes)){
    if(any(grepl(galliformes[[i]], x))){
      x <- names(galliformes)[[i]]
    }
  }
  
  return(x)
}

FormatStrigiformes <- function(x){
  strigiformes <- c(
    "owl sp." = "strigiformes \\(owl sp\\.\\)|owl sp\\.|^owl$",
    "short-eared owl" = "asio flammeus|short[- ]ear(ed){0,1} owl",
    "eurasian eagle-owl" = "bubo bubo|^(eurasian ){0,1}eagle[- ]owl$", 
    "spotted eagle-owl" = "bubo africanus|^spotted eagle[- ]owl$", 
    "tawny owl" = "strix aluco|tawny owl|towny owel",
    "long-eared owl" = "asio otus|long[- ]ear(ed){0,1} owl",
    "little owl" = "athene noctua|little owl",
    "barn owl" = "tyto alba|barn owl",
    "great horned owl" = "bubo virginianus|great horn(ed){0,1} owl",
    "ural owl" = "strix uralensis|ural owl",
    "eurasian scops-owl" = "otus scops|eurasian scops[- ]owl",
    "snowy owl" = "bubo scandiacus|snowy owl",
    "western screech-owl" = "megascops kennicottii|western screech[- ]owl",
    "northern hawk owl" =  "surnia ulula|northern[- ]hawk[- ]owl",
    "great grey owl"=	'great[- ]gr[ea]y[- ]owl|strix nebulosa'
    
  )
  
  for (i in 1:length(strigiformes)){
    if(any(grepl(strigiformes[[i]], x))){
      x <- names(strigiformes)[[i]]
    }
  }
  
  return(x)
}


FormatPasseriformes <- function(x){
  passeriformes <-  c(
    "eurasian jackdaw" = "corvus monedula|eurasian jackdaw|western jackdaw",
    "eurasian jay" = "garrulus glandarius|eurasian jay",
    "eurasian magpie" = "pica pica|eurasian magpie", # magpies divided by region
    "house sparrow" = "passer domesticus|house sparrow",
    "crow/raven sp." = "corvidae \\(crow/raven sp\\.\\)|crow/raven sp\\.|jungle crow|^crow$|^raven$|american raven",
    "large-billed crow" = "corvus macrorhynchos|large[- ]bill(ed){0,1} crow",
    "eurasian blackbird" =	"turdus merula|(eurasian|european) blackbird",
    'eurasian nuthatch' = '(eurasian|european) nuthatch|sitta europaea',
    "song thrush" = "turdus philomelos|song thrush",
    "black-billed magpie" = "pica hudsonia|black-billed magpie", # magpies divided by region
    "american crow" = "corvus brachyrhynchos|(american|northwestern) crow",
    "acrocephalus sp." = "acrocephalus \\(acrocephalus sp\\.\\)|acrocephalus sp\\.|^reed warbler$",
    "crow sp." = "corvus \\(crow sp\\.\\)|crow sp\\.|corvus corvus",
    "great-tailed grackle" = "quiscalus mexicanus|great[- ]tail(ed){0,1} grackle",
    "eurasian bullfinch" = "pyrrhula pyrrhula|eurasian bullfinch",
    "grey-headed bullfinch" = "pyrrhula erythaca|grey-headed bullfinch",
    "redwing" = "turdus iliacus|redwing",
    "fieldfare" = "turdus pilaris|fieldfare",
    "bali myna" = "leucopsar rothschildi|bali (starling|myna)",
    "eurasian blue tit" = "cyanistes caeruleus|eurasian blue tit",
    "common chaffinch" = "fringilla coelebs|common chaffinch",
    "pin-tailed parrotfinch" = "erythrura prasina|pin-tailed parrotfinch",
    "european robin" = "erithacus rubecula|european robin",
    "great tit"	 = "parus major|great tit",
    "house sparrow"	 = "passer domesticus|house sparrow",
    "white wagtail" = 'motacilla alba|white wagtail',
    "spotted flycatcher" = 'muscicapa striata|spotted flycatcher',
    "bearded reedling" = "panurus biarmicus|bearded reedling",
    "common redstart" = "phoenicurus phoenicurus|common redstart",
    "island canary" = 'serinus canaria|island canary',
    "european starling" = 'sturnus vulgaris|european starling',
    "eurasian tree sparrow"	= "passer montanus|eurasian tree sparrow",
    "common hill myna" = "common hill myna|religiosa",
    'thrush sp.' = '(turdidae|thrush)( sp\\.){0,1}',
    "pycnonotidae sp." = 'bul {0,1}bul'
  )
  
  for (i in 1:length(passeriformes)){
    if(any(grepl(passeriformes[[i]], x))){
      x <- names(passeriformes)[[i]]
    }
  }
  
  return(x)
}





FormatSuliformes <- function(x){
  suliformes <- c(
    "cormorant sp." = "phalacrocoracidae \\(cormorant sp\\.\\)|cormorant sp\\.|^cormorant$",
    "great cormorant" = "phalacrocorax carbo|(great[- ]|white[- ]breast(ed){0,1}[- ])cormorant",
    "northern gannet" = "morus bassanus|nor{0,1}thern gannet", # gannets divided by region
    "cape cormorant" = "phalacrocorax capensis|cape cormorant",
    "sulid sp." = "sulidae \\(sulid sp\\.\\)|sulid sp\\.|^gannet$",
    "cape gannet" = "morus capensis|cape gannet",
    "peruvian booby" = "sula variegata|peruvian booby",
    "guanay cormorant" = "leucocarbo bougainvillii|guanay cormorant",
    "neotropic cormorant" = "phalacrocorax brasilianus|neotropic cormorant",
    "double-crested cormorant" = "phalacrocorax auritus|double[- ]crest(ed){0,1} cormorant",
    "brown booby" = "sula leucogaster|brown booby",
    'magnificent frigatebird' = 'magnificent frigatebird|fregata[ -]magnificens'
  )
  
  for (i in 1:length(suliformes)){
    if(any(grepl(suliformes[[i]], x))){
      x <- names(suliformes)[[i]]
    }
  }
  
  return(x)
}


FormatGruiformes <- function(x){
  gruiformes <-  c(
    "crane sp." = "grus sp\\.|crane sp\\.|^crane$|black-throated crane",
    "great egret" = "ardea alba|great egret",
    "hooded crane" = "grus monacha|hooded crane",
    "blue crane" = "grus paradisea|blue crane",
    "coot sp." = "fulica sp\\.|coot sp\\.|^coot$",
    "rail/crake sp." = "rallidae sp\\.|rail/crake sp\\.|^moorhen$",
    "red-crowned crane" = "grus japonensis|red[- ]crown(ed){0,1}[- ]crane",
    "white-naped crane" = "grus vipio|white[- ]nap(ed){0,1}[- ]crane",
    'eurasian coot' = '(common|eurai{0,1}si{0,1}an|australian) coot|fulica atra'
    )

  
  for (i in 1:length(gruiformes)){
    if(any(grepl(gruiformes[[i]], x))){
      x <- names(gruiformes)[[i]]
    }
  }
  
  return(x)
}


FormatPodicipediformes <- function(x){
  podicipediformes <- c("grebe sp." = "podicipedidae sp\\.|grebe sp\\.|^grebe$",
                        "great crested grebe" = "podiceps cristatus|((great ){0,1}crested|(g c)) grebe",
                        "little grebe" = "tachybaptus ruficollis|little grebe",
                        "eared grebe" = "podiceps nigricollis|eared grebe")
  
  for (i in 1:length(podicipediformes)){
    if(any(grepl(podicipediformes[[i]], x))){
      x <- names(podicipediformes)[[i]]
    }
  }
  
  return(x)
}


FormatFalconiformes <- function(x){
  falconiformes <- c("peregrine falcon" = "falco peregrinus|peregrine falcon",
                     "falcon sp." = "falco sp\\.|falcon sp\\.|^falcon$",
                     'eurasian kestrel' = '(common|european|eurasian|old[- ]world) kestr{0,1}el|falco tinnunculus',
                     'american kestrel' = 'american kestrel|falco sparverius'
                     )
  
  
  for (i in 1:length(falconiformes)){
    if(any(grepl(falconiformes[[i]], x))){
      x <- names(falconiformes)[[i]]
    }
  }
  
  return(x)
}


FormatCiconiiformes <- function(x){
  ciconiiformes <- c(
    "white stork" = "ciconia ciconia|european white stork",
    "stork sp." = "ciconia sp\\.|stork sp\\.|^stork$"
  )
  
  for (i in 1:length(ciconiiformes)){
    if(any(grepl(ciconiiformes[[i]], x))){
      x <- names(ciconiiformes)[[i]]
    }
  }
  
  return(x)
}


FormatRheiformes <- function(x){
  rheiformes <- c(
    "greater rhea" = "rhea americana|greater rhea|^rhea$"
  )
  
  for (i in 1:length(rheiformes)){
    if(any(grepl(rheiformes[[i]], x))){
      x <- names(rheiformes)[[i]]
    }
  }
  
  return(x)
}


FormatColumbiformes <- function(x){
  columbiformes <- c(
    "common wood-pigeon" = "columba palumbus|common wood[-] pigeon",
    "rock pigeon" = 'rock pigeon|columba livia',
    "eurasian collared-dove" = "eurasian collared[- ]dove|streptopelia decaocto",
    "pigeon/dove sp." = "columbidae sp\\.|pigeon/dove sp\\.|^pigeon$|^dove$|columbidae|columba sp\\."
  )
  
  for (i in 1:length(columbiformes)){
    if(any(grepl(columbiformes[[i]], x))){
      x <- names(columbiformes)[[i]]
    }
  }
  
  return(x)
}


FormatPsittaciformes <- function(x){
  psittaciformes <-  c(
    "amazona sp." = "amazona sp\\.|amazona sp\\.|amazon parrot",
    "large macaw sp." = "large macaw sp\\.|large macaw sp\\.|catalina macaw",
    "parrot sp." =	"psittaciformes sp\\.|parrot sp\\.|^parrot$|ornamental bird"
  )
  
  for (i in 1:length(psittaciformes)){
    if(any(grepl(psittaciformes[[i]], x))){
      x <- names(psittaciformes)[[i]]
    }
  }
  
  return(x)
}


FormatCathartiformes <- function(x){
  cathartiformes <-c(
    "black vulture" = "coragyps atratus|black vulture"

  )
  
  for (i in 1:length(cathartiformes)){
    if(any(grepl(cathartiformes[[i]], x))){
      x <- names(cathartiformes)[[i]]
    }
  }
  
  return(x)
}


FormatProcellariiformes <- function(x){
  procellariiformes <- c(
    "northern giant-petrel" = "macronectes halli|northern giant[- ]petrel"
  )
  
  for (i in 1:length(procellariiformes)){
    if(any(grepl(procellariiformes[[i]], x))){
      x <- names(procellariiformes)[[i]]
    }
  }
  
  return(x)
}


FormatSphenisciformes <- function(x){
  sphenisciformes <- c("african penguin" = "spheniscus demersus|african penguin",
                       "humboldt penguin" = "spheniscus humboldti|humboldt penguin")
  
  for (i in 1:length(sphenisciformes)){
    if(any(grepl(sphenisciformes[[i]], x))){
      x <- names(sphenisciformes)[[i]]
    }
  }
  
  return(x)
}


FormatCasuariiformes <- function(x){
  casuariiformes <- c( "emu" = "dromaius novaehollandiae|^emu$")
  
  for (i in 1:length(casuariiformes)){
    if(any(grepl(casuariiformes[[i]], x))){
      x <- names(casuariiformes)[[i]]
    }
  }
  
  return(x)
}


FormatStruthioniformes <- function(x){
  struthioniformes <- c("common ostrich" = "struthio camelus|^(common ){0,1}ostrich$")
  
  for (i in 1:length(struthioniformes)){
    if(any(grepl(struthioniformes[[i]], x))){
      x <- names(struthioniformes)[[i]]
    }
  }
  
  return(x)
}


FormatPhoenicopteriformes <- function(x){
  phoenicopteriformes <- c("flamingo sp." = "phoenicopterus sp\\.|flamingo sp\\.|^flamingo$")
  
  for (i in 1:length(phoenicopteriformes)){
    if(any(grepl(phoenicopteriformes[[i]], x))){
      x <- names(phoenicopteriformes)[[i]]
    }
  }
  
  return(x)
}


FormatCaprimulgiformes <- function(x){
  caprimulgiformes <- c("eurasian nightjar" =	"caprimulgus europaeus|eurasian nightjar")
  
  for (i in 1:length(caprimulgiformes)){
    if(any(grepl(caprimulgiformes[[i]], x))){
      x <- names(caprimulgiformes)[[i]]
    }
  }
  
  return(x)
}
###################################################################################################


FormatBird<- function(x){
  x <- FormatAccipitriformes(x) 
  x <- FormatAnseriformes(x)
  x <- FormatCasuariiformes(x) 
  x <- FormatCathartiformes(x) 
  x <- FormatCharadriiformes(x) 
  x <- FormatCiconiiformes(x) 
  x <- FormatColumbiformes (x) 
  x <- FormatFalconiformes(x) 
  x <- FormatGalliformes(x) 
  x <- FormatGruiformes(x) 
  x <- FormatPelecaniformes(x) 
  x <- FormatPasseriformes(x) 
  x <- FormatPhoenicopteriformes(x)
  x <- FormatPodicipediformes(x) 
  x <- FormatProcellariiformes(x) 
  x <- FormatPsittaciformes(x) 
  x <- FormatRheiformes(x) 
  x <- FormatSphenisciformes(x) 
  x <- FormatStrigiformes(x) 
  x <- FormatStruthioniformes(x) 
  x <- FormatSuliformes(x) 
  x <- FormatCaprimulgiformes(x)
  x <- gsub('^birds{0,1}$|avian|.*wild {0,1}birds{0,1}.*|ac{0,1}quatic birds{0,1}', 'aves', x)
  x <- gsub("^env{0,1}(iron{0,1}ment{0,1}){0,1}.*|^water$", 'environment', x)

  return(x)  
}


###################################################################################################
# Unused bird orders

#FormatCoraciformes
#FormatBucerotiformes
#FormatTrogoniformes
#FormatLeptosomiformes
#FormatColilformes
#FormatOpisthocomiformes
#FormatGaviiformes
#FormatPterocliformes
#FormatCucliformes
#FormatOtidiformes
#FormatMusophagiformes
#FormatApodiformes
#FormatCaprimulgiformes



#FormatApterygiformes
#FormatTinamiformes