#ifndef __DD_H__
#define __DD_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector.h>
#include "mrf.h"
#include "LinkedBlockList.h"

void print2DVect(vector<vector<float> > &vectorToPrint, int size1, int size2);
void print2DVect(vector<vector<int> > &vectorToPrint, int size1, int size2);

void null2DVect(vector<vector<float> > &vectToNullify, int size1, int size2);

void null3DVect(vector<vector<vector<float> > > &vectToNullify, int size1, int size2, int size3);

void copyVect(vector<MRF::Label> &toCopy, MRF::Label* &toFill);

class DD : public MRF{
public:
    EnergyVal smoothnessEnergy();
    EnergyVal dataEnergy();
    MRF::EnergyVal dualEnergy();
    void writeEnergies();
    
    bool strongAgreement;
    
    DD(int width, int height, int nLabels, EnergyFunction *eng);
    DD(int nPixels, int nLabels,EnergyFunction *eng);
    ~DD();
    
    void setNeighbors(int pix1, int pix2, CostVal weight);
    Label getLabel(int pixel){return(m_answer[pixel]);};
    void setLabel(int pixel,Label label){m_answer[pixel] = label;};
    Label* getAnswerPtr(){return(m_answer);};
    void setParameters(int /*numParam*/, void * /*param*/){printf("No optional parameters to set"); exit(1);}
    
    void minSum(vector<vector<float> > &chain, vector<Label> &answer, int indexChain);
    
    bool doTreesAgree();
    void printHorizontalAnswers();
    void printVerticalAnswers();
    
    
protected:
    void initializeAlg();
    void clearAnswer();
    void optimizeAlg(int nIterations);
    
    void setData(DataCostFn dcost); 
    void setData(CostVal* data);    
    void setSmoothness(SmoothCostGeneralFn cost);
    void setSmoothness(CostVal* V);
    void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    void setCues(CostVal* hCue, CostVal* vCue); 
    
private:
    Label *m_answer;
    vector<Label> m_answer_horizontal;
    vector<Label> m_answer_vertical;
    CostVal *m_V;
    CostVal *m_D;
    
    CostVal *m_horizWeights;
    CostVal *m_vertWeights;
    
    int m_nChains;
    
    vector<vector<vector<float> > > dualHorizontalChains;
    vector<vector<vector<float> > > dualVerticalChains;
    
    
    DataCostFn m_dataFn;
    SmoothCostGeneralFn m_smoothFn;
    
    typedef struct NeighborStruct{
        int     to_node;
        CostVal weight;
    } Neighbor;
    
    LinkedBlockList *m_neighbors;
};

#endif