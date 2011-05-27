      integer NDAT1,NDAT2,NDAT3,NDAT4,NDAT5,NDAT6,NDAT7,NDAT8,NDAT9
      integer NDAT10,NDAT11,NDAT12,NDAT13,NDAT14,NDAT15,NDAT16,NDAT17
      integer NDAT18
      integer NDAT19,NDAT20
      integer NDAT21,NDAT22
      integer NDAT23,NDAT31
      integer NDAT24   ! added for H1 Low Ep dataset 05/08/2010
      integer NDAT25  ! ZEUS prel
      integer NDAT26, NDAT27 ! single diff H1	
      integer NDAT32,NDAT33,NDAT34,NDAT35,NDAT36,NDAT37,NDAT38,NDAT39
      integer NDAT40,NDAT41,NDAT42
cv...
      integer NDAT61
      integer NDAT62,NDAT63
      integer NDAT64      

cv+sh 
      integer NDAT65
      integer NDAT66

      integer NDAT67
      integer NDAT68
      integer NDAT69
      integer NDAT70

cv fixed target
	integer NDAT91
	integer NDAT92
	integer NDAT93
	integer NDAT94
	integer NDAT95
	integer NDAT96

cv atlas
        integer NDAT100
        integer NDAT101
        integer NDAT102
        integer NDAT103
        integer NDAT104
        
cv tevtraon
	integer NDAT111
        integer NDAT112

          
      integer NTOT

      parameter (NDAT1 = 45)
      parameter (NDAT2 = 80)
      parameter (NDAT3 = 130)
      parameter (NDAT4 = 25)
      parameter (NDAT5 = 126)
      parameter (NDAT6 = 28)
      parameter (NDAT7 = 13)
      parameter (NDAT8 = 147)
      parameter (NDAT9 = 28)
c bulk00 
      parameter (NDAT10 = 99)
c low q2 comb 95-00
      parameter (NDAT11 = 149)
      parameter (NDAT12 = 262)

          
      parameter (NDAT13 = 147)
      parameter (NDAT14 = 150)
      parameter (NDAT15 = 29)
      parameter (NDAT16 = 29)

      parameter (NDAT17 = 21)
      parameter (NDAT18 = 26)
      parameter (NDAT19 = 27)
      parameter (NDAT20 = 26)
C -- F2cc and F2bb, low Q2
      parameter (NDAT21 = 10)
      parameter (NDAT22 = 10)
C -- F2cc D* data, low Q2
      parameter (NDAT23 = 12)
      parameter (NDAT24 = 136)   ! added 05/08/2010; changed 13/08/2010
      parameter (NDAT25 = 90)
      parameter (NDAT26 = 17)
      parameter (NDAT27 = 17)

C------------------------------------------          
C -- Jet data
      parameter (NDAT31 = 16)
      parameter (NDAT32 = 19)
      parameter (NDAT33 = 16)
      parameter (NDAT34 = 24)
      parameter (NDAT35 = 24)
      parameter (NDAT36 = 24)
      parameter (NDAT37 = 28)
      parameter (NDAT38 = 30)
      parameter (NDAT39 = 6)

      parameter (NDAT40 = 70)
      parameter (NDAT41 = 30)
      parameter (NDAT42 = 30)

      


C------------------------------------------   
cv .. combo data
      parameter (NDAT61 = 145)
      parameter (NDAT62 = 383)
      parameter (NDAT63 = 34)
      parameter (NDAT64 = 34)
      parameter (NDAT65 = 41)
      parameter (NDAT66 = 230)
      parameter (NDAT67 = 182)
      parameter (NDAT68 = 561)
      parameter (NDAT69 = 41)
      parameter (NDAT70 = 39)
c      parameter (NDAT91 = 661)        
      parameter (NDAT91 = 211)    
      parameter (NDAT92 = 158)    
      parameter (NDAT93 = 178)    


cv atlas


      parameter (NDAT100 = 7)    
      parameter (NDAT101 = 11)
      parameter (NDAT102 = 11)
      parameter (NDAT103 = 11)
      parameter (NDAT104 = 11)    

cv tevatron          
      parameter (NDAT111 = 12)        
      parameter (NDAT112 = 13)        
          
      parameter(NTOT = NDAT1+NDAT2+NDAT3+NDAT4+NDAT5+NDAT6+
     +                 NDAT7+NDAT8+NDAT9+NDAT10+NDAT11+NDAT12+
     +                 NDAT13+NDAT14+NDAT15+NDAT16+NDAT17+NDAT66+
     +                 NDAT18+NDAT19+NDAT20+NDAT21+NDAT22+NDAT23+
     +                 NDAT25+NDAT26+NDAT27+
     +                 NDAT31+NDAT32+NDAT33+NDAT34+NDAT35+NDAT36+
     +                 NDAT37+NDAT38+NDAT39+
     +                 NDAT40+NDAT41+NDAT42+
     +                 NDAT61+NDAT62+NDAT63+NDAT64+NDAT65+
     +                 NDAT68+NDAT69+NDAT70+NDAT67+NDAT91+NDAT92+
     +                 NDAT93+NDAT100+NDAT101+NDAT102+NDAT103+NDAT104+           
     +                 NDAT111+NDAT112)
