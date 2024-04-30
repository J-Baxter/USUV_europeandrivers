FormatMammal <- function(x){
  mammals <- c(
    "red fox" = "(red|blue|silver) fox|^fox$|vulpes vulpes",
    "brown rat" = "r{0,1}attus norvegicus",
    'harbour seal'= "harbou{0,1}r seal",
    'ferret'= "mustela furo|^ferret$",
    'european polecat'= 'mustela putorius|european polecat',
    'beech marten'= "(stone|beech) marten",
    "fisher"= "pekania pennanti|^fisher$",
    "arctic fox"= "arctic[ -]fox|vulpes lagopus",
    'human' = '^humans{0,1}$|homo sapiens',
    'lion' = '^lion$|leo panthera',
    'eared seal sp.' = '^sea lion$',
    'american black bear' = '^black bear$|ursus americanus',
    'true seal sp.' = '^seal$|seal sp\\.',
    'atlantic grey seal' = 'gr[ea]y seal',
    'bear sp.'= "^bear$|bear sp\\.",
    'mustelid sp.'= '^mink$|wild mink|mink sp\\.|polecat|^otter$|badger',
    'domestic ferret' = '^ferret$|mustela furo', 
    'common otter' = 'lutra {0,1}lutra|common otter',
    "dolphin sp."= "^dolphin$",
    'common dolphin' = '(short[- ]beaked|long[- ]beaked) common dolphin|^common dolphin$',
    "porpoise sp."= "^porpoise$",
    "lynx"= "lynx sp\\.|^lynx$",
    'skunk sp.'= '^skunk$|skunk sp\\.',
    'feline sp.'= '^cat$|domestic cat|feline',
    'canid sp.'= 'canine',
    "shrew sp." = '(soricidae|shrew)( sp\\.){0,1}',
    "coucha rat" = "mastomys natalensis|coucha rat",
    "black rat"="rattus rattus|black rat",
    'common pipistrelle' = 'pipistrellus pipistrellus|common pipistrelle'
    )
  
  for (i in 1:length(mammals)){
    if(any(grepl(mammals[[i]], x))){
      x <- names(mammals)[[i]]
    }
  }
  
  return(x)
}

