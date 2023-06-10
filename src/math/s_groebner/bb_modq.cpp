
#include <s_groebner/bb_modq.h>
#include <val_utils.h>
#include <MyTime.h>


namespace bb_modq
{

int nGrest=0,nGrestold=0,nGall=0;



void update(val::s_polynom<val::modq> &f,val::Glist< val::s_polynom<val::modq> > &G,
            val::Glist< val::s_polynom<val::modq> > &Grest,val::Glist< common_bb::spair > &lspair,int &m)
{
 int i,j;
 val::s_expo t1;
 int *h_done;
 val::s_polynom<val::modq> *h;
 val::GlistManipulator<val::s_polynom<val::modq> > ItG,ItG2;
 val::GlistIterator<common_bb::spair > ItPair;

 if (f.iszero()) return;
 if (G.isempty()) {
	 G.sinsert(std::move(f));
	 m=1;
	 return;
 }

 h_done = new int[m];
 for (i=0;i<m;i++) h_done[i]=0;

 for (i=0,ItG.settohead(G);i<m;i++,ItG.moveactual()) { //  1.while-loop
	 if (h_done[i]) continue;
	 if (common_bb::lcmdis(f.LT(),ItG.getelement().LT(),t1)) {// keeping still disjoint pairs (h,P[i])
		 h_done[i]=2;
		 continue;
	 }
	 for (j=0,ItG2.settohead(G);j<m;j++,ItG2.moveactual()) {
		 if (j==i || h_done[j]==1) continue;
		 if (val::lcm(f.LT(),ItG2.getelement().LT())|t1) {
			 h_done[i]=1;
			 break;
		 }
	 }
 }

 for (lspair.resetactual();lspair.actualvalid();lspair.moveactual()) {
	 if (lspair.getelement().s_done) continue;
	 t1 = val::lcm(lspair.getelement().fmodq->LT(),lspair.getelement().gmodq->LT());
	 if ( (f.LT()|t1) && (val::lcm(lspair.getelement().fmodq->LT(),f.LT())!=t1)
		 && (val::lcm(f.LT(),lspair.getelement().gmodq->LT())!=t1) ) {
		 lspair.actualvalue().s_done=1;
	 }
 }

 // Ordered insertion in G: G.head!=NULL.
 G.resetactual();
 if (f.LT()< G.actualvalue().LT()) {
     G.inserttohead(std::move(f));
     G.resetactual();
     h=&G.actualvalue();
     j=0;
 }
 else {
	 for (G.resetactual(),j=1;G.nextvalid() && G.getnextelement().LT() < f.LT();G.moveactual(),j++);
	 G.insertnexttoactual(std::move(f));
	 G.moveactual();
	 h=&G.actualvalue();
 }
 m++;

 // Insertion in dPairs:
 int k=0;
 for (i=0,ItG.settohead(G);i<m;i++,ItG.moveactual(),k++) {
	 if (i==j) {k--;continue;}
	 if (h_done[k]) continue;
	 lspair.sinsert(common_bb::spair(&ItG.actualvalue(),h));
 }

 while (G.nextvalid()) {
	 if (h->LT() | G.getnextelement().LT()) {
         G.movenextelementtohead(Grest);
         m--;
         nGrest++;
	 }
	 else G.moveactual();
 }
 G.resetactual();

  // Eventually release memory:
  if (nGrest-nGrestold>=100) {
	  int is;//nelm=0;
	  for (ItG.settohead(Grest);ItG;ItG++) {
		  if (ItG().iszero()) continue;
		  is=1;
		  for (ItPair.settohead(lspair);ItPair;ItPair++) {
			  if (ItPair().s_done) continue;
			  if (ItPair().fmodq == &ItG() || ItPair().gmodq == &ItG()) {
				  is=0;
				  break;
			  }
		  }
		  if (is) {
			  nGrest--;
			  //nelm++;
			  ItG().del();
		  }
	  }
	  nGrestold=nGrest;

	  // Delete unnecessary pairs in lspair:
	  lspair.resetactual();
	  while (!lspair.isempty() && (lspair.actualvalue().s_done)) {
		 lspair.skiphead();
	  }

	  if (!lspair.isempty()) {
		  for (lspair.resetactual();lspair.nextvalid();) {
			  if (lspair.getnextelement().s_done) {
				  lspair.skipnextelement();
			  }
			   else lspair.moveactual();
		  }
	  }
	  lspair.resetlast();
	  lspair.resetactual();
	  //
  }
  delete[] h_done;
}



int Groebner(val::Glist< val::s_polynom<val::modq> > &G,int comment)
{
 int m=0,anzspair=0;
 val::s_polynom<val::modq> h;
 val::Glist< val::s_polynom<val::modq> > H,Grest;
 val::Glist< common_bb::spair > lspair;
 std::string s;

 nGrest=0,nGrestold=0,nGall=0;
 val::s_polynom<val::modq>::nreduction=0;

 // Set G and H:
 if (G.isempty()) return 0;
 H.copyanddelete(G);


 while ((!H.isempty()) || (!lspair.isempty())) {
	 //Choose minimal element:
	 if (!H.isempty() && !lspair.isempty()) {
         H.resetactual();lspair.resetactual();
		 if (val::lcm(lspair.getelement().fmodq->LT(),lspair.getelement().gmodq->LT())< (H.getelement().LT())) {
			 anzspair++;
			 h=common_bb::spol(*(lspair.getelement().fmodq),*(lspair.getelement().gmodq));
			 lspair.skiphead();
		 }
		 else {
			 h=std::move(H.actualvalue());
			 H.skiphead();
         }
	 }
	 else if (!H.isempty()) { // lspair empty
		 H.resetactual();
		 h=std::move(H.actualvalue());
		 H.skiphead();
	 }
	 else { // H empty
		 lspair.resetactual();
		 anzspair++;
         h=common_bb::spol(*(lspair.getelement().fmodq),*(lspair.getelement().gmodq));
		 lspair.skiphead();
	 }
	 h.reduction(G,0);
	 if (!h.iszero()) {
		 update(h,G,Grest,lspair,m);
		 nGall++;
	 }

	 lspair.resetactual();
	 while(!lspair.isempty() && lspair.getelement().s_done) lspair.skiphead();
 }
 Grest.dellist();
 common_bb::interredBasis(G);
 if (comment) {
    s+="\nNumber of used critical pairs: "+val::ToString(anzspair);
    s+="\nTotal number of all computed polynomials: " + val::ToString(nGall);
    common_bb::WriteText(s);
 }
 return m;
}


int bbgccmain(const std::string &argv,val::Glist< val::s_polynom<val::modq> > &G)
{
 int m,zeit;
 std::string s;

 if (!common_bb::readfromfile(argv,G)) return 0;

 val::ChronoClass Chrono;
 m=Groebner(G);
 zeit=Chrono();
 common_bb::Clear();
 s= "\n\nG has " + val::ToString(m) + " elements."
  + "\nNumber of monomials: " + val::ToString(val::s_polynom<val::modq>::getmnumber())
  + "\nReductions: " + val::ToString(val::s_polynom<val::modq>::nreduction)
  + "\nTime in sec.: " + val::ToString(zeit);

 common_bb::WriteText(s);
 nGrest=0;nGrestold=0;nGall=0;
 val::s_polynom<val::modq>::nreduction=0;

 return m;
}


} // end namespace bb_modq
