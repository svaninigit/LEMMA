
   // Declaration of leaf types
   Int_t           iev;
   Int_t           nhits;
   // MC variables as defined by M. Dreucci 
   Int_t           subdet[100];   //[nhits]
   Int_t           idp[100];   //[nhits]
   Int_t           ipar[100];   //[nhits]
   Int_t           itrack[100];   //[nhits]
   Double_t        time[100];   //[nhits]
   Double_t        xh[100];   //[nhits]
   Double_t        yh[100];   //[nhits]
   Double_t        zh[100];   //[nhits]
   Double_t        p[100];   //[nhits]
   Double_t        pxh[100];   //[nhits]
   Double_t        pyh[100];   //[nhits]
   Double_t        pzh[100];   //[nhits]
   Double_t        xv[100];   //[nhits]
   Double_t        yv[100];   //[nhits]
   Double_t        zv[100];   //[nhits]
   Double_t        kinev[100];   //[nhits]
   Double_t        pxvdir[100];   //[nhits]
   Double_t        pyvdir[100];   //[nhits]
   Double_t        pzvdir[100];   //[nhits]
   Int_t           pro[100];   //[nhits]
   Int_t           istep[100];   //[nhits]
   Int_t           inextstep[100];   //[nhits]
   // Additional variables (Tue May  9 2017)
   // Si cluster lenght
   Int_t           clulen[100];   //[nhits]
   // Gamma calorimeter (X=Y=0)
   Double_t        gcal[64];
   // Positron calorimeter (at X<0)
   Double_t        pcal[64];
   // Additional calorimeter 
   Double_t        calo[64];
