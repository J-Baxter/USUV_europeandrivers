FormatTicks <- function(x){
  ticks <- c(
    # hard ticks
    'ixodidae spp.' = 'hard tick',
   'ixodes ricinus' = 'ixodes ricinus',
   'ixodes persulcatus' = 'ixodes persulcatus',
   'ixodes frontalis' = 'ixodes frontalis',
   'hyalomma lusitanicum' = 'hyalomma lusitanicum',
   'hyalomma marginatum' = 'hyalomma marginatum',
   'dermacentor reticulatus' = 'dermacentor reticulatus',
   'rhipicephalus sanguineus' = 'rhipicephalus sanguineus',
   'ixodes spp.' = 'ixodes4',
   'dermacentor spp.' = 'dermacentor$',
   'hyalomma spp.' = 'hyalomma$',
   'rhipicephalus spp.' = 'rhipicephalus$',
  
   
   # soft ticks
   'argasidae spp.' = 'soft tick',
   'carios erraticus' = '(carios|ornithodours) erraticus',
   'ornithodours' = '(ornithodours|carios)$',
   
   'ixodida spp.' = '^tick' )
  
  for (i in 1:length(ticks)){
    if(any(grepl(ticks[[i]], x))){
      x <- names(ticks)[[i]]
    }
  }
  
  return(x)
}

