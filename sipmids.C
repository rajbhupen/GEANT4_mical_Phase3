void sipmids(){



      unsigned int id=0;
      int cmax,bmax;
      
      for(int ia=1; ia<8; ia++){//locno
	id++;
	id<<=2;
				
	//	cmax = (ia==1) ? 88:40;
	bmax = (ia==1) ? 4 : 3;
	if(ia==1){cmax=88;}
	else if(ia>5){cmax = 8;}
	else{cmax = 40;}
	for(int jb =0; jb<bmax; jb++) {//layerno
	  if(jb==0){ id+=0;} else{id++;}
	  id<<=7;
	  
	  for(int kc =0; kc<cmax; kc++){//stripno
	    if(kc==0){id+=0;} else{	id++; }
	    id<<=2;
	    
	    for(int med=0; med<4; med++){//sipmno
	      if(med==0){id+=0;} else{id++;	} 
	      
	      cout<<"id " <<id<<"  "<<ia<<"  "<<jb<<"  "<<kc<<"  " <<med<<endl;
	   

	    }//sipm no
	    id-=3;
	    id>>=2;
	  }//stripno
	  id-=(cmax-1);
	  id>>=7;
	  
	  //cout<<endl;
	}//layerno
	id-=(bmax-1);
	id>>=2;
	//cout<<endl<<endl;
      }//locno
      






}
