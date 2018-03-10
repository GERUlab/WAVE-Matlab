
((WC(k)-WCB(k))./dt+rtex(1:length(WC)))=(-wat_flxs(k)...
      +wat_flxs(k+1))./dx
  
      
  wat_flxs(k+1)= ((WC(k)-WCB(k))./dt+rtex(1:length(WC)))*dx +wat_flxs(k)