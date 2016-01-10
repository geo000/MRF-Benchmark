#include <string.h>
#include <assert.h>
#include <limits.h>
#include <iostream> 
#include <fstream>
#include "mrf.h"
#include "DD.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const int Inf = numeric_limits<int>::max();

#define pixel(x,y)  x+y*m_width
#define index(x,y,label)  pixel(x,y)*m_nLabels + label
#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]

DD::DD(int width, int height, int nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
}
DD::DD(int nPixels, int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng)
{
}

DD::~DD()
{ 
    delete[] m_answer;
    if (!m_grid_graph) delete[] m_neighbors;
}

void DD::initializeAlg() {
    int pix(0);
    // initialize answer
    m_answer = (Label *) new Label[m_nPixels];
    
    // initialize constants
    strongAgreement = false;
    m_nChains = m_width + m_height;
    
    // initialize chains
    m_answer_horizontal.resize(m_nPixels);
    m_answer_vertical.resize(m_nPixels);
    null3DVect(dualHorizontalChains, m_height, m_width, m_nLabels);
    null3DVect(dualVerticalChains, m_width, m_height, m_nLabels);
    
    for (int iCol=0; iCol<m_height; iCol++ ) {
        for (int iLine=0; iLine<m_width; iLine++) {
            for (Label iLabel=0; iLabel<m_nLabels; iLabel++) {
                pix = iCol*m_width + iLine;
                dualHorizontalChains[iCol][iLine][iLabel] = m_D(pix,iLabel);
            }
        }
    }
    
    for (int iLine=0; iLine<m_width; iLine++) {
        for (int iCol=0; iCol<m_height; iCol++ ) {
            for (Label iLabel=0; iLabel<m_nLabels; iLabel++) {
                pix = iCol*m_width + iLine;
                dualVerticalChains[iLine][iCol][iLabel] = m_D(pix,iLabel);
            }
        }
    }
}


void DD::clearAnswer()
{
    memset(m_answer, 0, m_nPixels*sizeof(Label));
}

void DD::setNeighbors(int pixel1, int pixel2, CostVal weight)
{
    assert(!m_grid_graph);
    assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);
    
    
    Neighbor *temp1 = (Neighbor *) new Neighbor;
    Neighbor *temp2 = (Neighbor *) new Neighbor;
    
    if ( !temp1 || ! temp2 ) {printf("\nNot enough memory, exiting");exit(0);}
    
    temp1->weight = weight;
    temp1->to_node = pixel2;
    
    temp2->weight  = weight;
    temp2->to_node = pixel1;
    
    m_neighbors[pixel1].addFront(temp1);
    m_neighbors[pixel2].addFront(temp2);
}

MRF::EnergyVal DD::dualEnergy() {
    EnergyVal eng = (EnergyVal) 0;
    
    // Unary energy
    for (int iLine = 0; iLine < m_width; iLine++ ) {
        for (int iCol=0; iCol < m_height; iCol++) {
            eng += dualVerticalChains[iLine][iCol][m_answer_vertical[pixel(iLine, iCol)]] + dualHorizontalChains[iCol][iLine][m_answer_horizontal[pixel(iLine, iCol)]];
        }
    }
    
    // Binary energy
    for (int iLine = 0; iLine < m_width-1; iLine++ ) {
        for (int iCol=0; iCol < m_height-1; iCol++) {
            eng += m_V(m_answer_horizontal[pixel(iLine, iCol)], m_answer_horizontal[pixel(iLine+1, iCol)]) + m_V(m_answer_vertical[pixel(iLine, iCol)], m_answer_vertical[pixel(iLine, iCol+1)]);
        }
    }
    
    return eng;
}

MRF::EnergyVal DD::smoothnessEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    EnergyVal weight;
    int x,y,pix,i;
    
    if ( m_grid_graph )
    {
        if ( m_smoothType != FUNCTION  )
        {
            for ( y = 0; y < m_height; y++ )
                for ( x = 1; x < m_width; x++ )
                {
                    pix    = x+y*m_width;
                    weight = m_varWeights ? m_horizWeights[pix-1] :  1;
                      eng = eng + m_V(m_answer[pix],m_answer[pix-1])*weight;
                }
            
            for ( y = 1; y < m_height; y++ )
               for ( x = 0; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    weight = m_varWeights ? m_vertWeights[pix-m_width] :  1;
                    eng = eng + m_V(m_answer[pix],m_answer[pix-m_width])*weight;
                }
        }
        else
        {
            for ( y = 0; y < m_height; y++ )
                for ( x = 1; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    eng = eng + m_smoothFn(pix,pix-1,m_answer[pix],m_answer[pix-1]);
                }
            
            for ( y = 1; y < m_height; y++ )
                for ( x = 0; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    eng = eng + m_smoothFn(pix,pix-m_width,m_answer[pix],m_answer[pix-m_width]);
                }
        }
    }
    else
    {
        
        Neighbor *temp; 
        
        if ( m_smoothType != FUNCTION  )
        {
            for ( i = 0; i < m_nPixels; i++ )
                if ( !m_neighbors[i].isEmpty() )
                {
                    m_neighbors[i].setCursorFront();
                    while ( m_neighbors[i].hasNext() )
                    {
                        temp = (Neighbor *) m_neighbors[i].next();
                        
                        if ( i < temp->to_node )
                            eng = eng + m_V(m_answer[i],m_answer[temp->to_node])*(temp->weight);
                    }
                }
        }
        else
        {
            for ( i = 0; i < m_nPixels; i++ )
                if ( !m_neighbors[i].isEmpty() )
                {
                    m_neighbors[i].setCursorFront();
                    while ( m_neighbors[i].hasNext() )
                    {
                        temp = (Neighbor *) m_neighbors[i].next();
                        if ( i < temp->to_node )
                            eng = eng + m_smoothFn(i,temp->to_node, m_answer[i],m_answer[temp->to_node]);
                    }
                }
        }
        
    }
    
    return(eng);
;}

MRF::EnergyVal DD::dataEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    
    
    if ( m_dataType == ARRAY) 
    {
        for ( int i = 0; i < m_nPixels; i++ )
            eng = eng + m_D(i,m_answer[i]);
    }
    else
    {
        for ( int i = 0; i < m_nPixels; i++ )
            eng = eng + m_dataFn(i,m_answer[i]);
    }
    return(eng);
    
}

void DD::writeEnergies() {
    ofstream fichier("energy.txt", ios::out | ios::app);
    fichier << smoothnessEnergy() + dataEnergy() << ", " << dualEnergy() << "\n";
    fichier.close();
}

void DD::setData(DataCostFn dcost)
{
    m_dataFn = dcost;
}

void DD::setData(CostVal* data)
{
    m_D = data;
}

void DD::setSmoothness(SmoothCostGeneralFn cost)
{
    m_smoothFn = cost;
}
void DD::setSmoothness(CostVal* V)
{
    m_V = V;
}

void DD::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    int i, j;
    CostVal cost;
    
    m_V = (CostVal *) new CostVal[m_nLabels*m_nLabels*sizeof(CostVal)];
    if (!m_V) { fprintf(stderr, "Not enough memory!\n"); exit(1); }
    
    
    for (i=0; i<m_nLabels; i++)
        for (j=i; j<m_nLabels; j++)
        {
            cost = (CostVal) ((smoothExp == 1) ? j - i : (j - i)*(j - i));
            if (cost > smoothMax) cost = smoothMax;
            m_V[i*m_nLabels + j] = m_V[j*m_nLabels + i] = cost*lambda;
        }
    
}

void DD::setCues(CostVal* hCue, CostVal* vCue)
{
    m_horizWeights = hCue;
    m_vertWeights  = vCue;
}


void DD::optimizeAlg(int nIterations)
{
    if (strongAgreement) {
        return;
    }
    
    float eta=1.0/(nIterations+1);
    
    // Compute optimal labeling
    for (int iChain=0; iChain<m_height; iChain++) {
        minSum(dualHorizontalChains[iChain], m_answer_horizontal, iChain);
    }
    for (int iChain=m_height; iChain<m_nChains; iChain++) {
        minSum(dualVerticalChains[iChain-m_height], m_answer_vertical, iChain);
    }
    
    // Update dual chains
    for (int iLine=0; iLine<m_width; iLine++) {
        for (int iCol=0; iCol<m_height; iCol++) {
            dualHorizontalChains[iCol][iLine][m_answer_horizontal[pixel(iLine, iCol)]] += eta*0.5;
            dualHorizontalChains[iCol][iLine][m_answer_vertical[pixel(iLine, iCol)]] -= eta*0.5;
            dualVerticalChains[iLine][iCol][m_answer_horizontal[pixel(iLine, iCol)]] -= eta*0.5;
            dualVerticalChains[iLine][iCol][m_answer_vertical[pixel(iLine, iCol)]] += eta*0.5;
        }
    }
    strongAgreement =  doTreesAgree();
    copyVect(m_answer_vertical, m_answer);
    writeEnergies();
}

void print2DVect(vector<vector<float> > &vectorToPrint, int size1, int size2) {
    for (int j=0; j<size2; j++) {
        for (int i=0; i<size1; i++) {
            cout << vectorToPrint[i][j] << " ";
        }
        cout << endl;
    }
}

void print2DVect(vector<vector<int> > &vectorToPrint, int size1, int size2) {
    for (int j=0; j<size2; j++) {
        for (int i=0; i<size1; i++) {
            cout << vectorToPrint[i][j] << " ";
        }
        cout << endl;
    }
}

void null2DVect(vector<vector<float> > &vectToNullify, int size1, int size2) {
    vectToNullify.resize(size1);
    for (int iSize1=0; iSize1 < size1; iSize1++) {
        vectToNullify[iSize1].resize(size2);
    }
}

void null3DVect(vector<vector<vector<float> > > &vectToNullify, int size1, int size2, int size3) {
    vectToNullify.resize(size1);
    for (int iSize1=0; iSize1<size1; iSize1++) {
        vectToNullify[iSize1].resize(size2);
        for (int iSize2=0; iSize2<size2; iSize2++) {
            vectToNullify[iSize1][iSize2].resize(size3);
        }
    }
}

void copyVect(vector<MRF::Label> &toCopy, MRF::Label* &toFill) {
    for (int iPix=0; iPix<toCopy.size(); iPix++) {
        toFill[iPix] = toCopy[iPix];
    }
}

bool DD::doTreesAgree() {
    for (int iPix=0; iPix<m_nPixels; iPix++) {
        if (m_answer_vertical[iPix] != m_answer_horizontal[iPix]) {
            return false;
        }
    }
    return true;
}

void DD::printHorizontalAnswers() {
    cout << "Printing hori :\n";
    for (int iCol=0; iCol<m_height; iCol++) {
        for (int iLine=0; iLine<m_width; iLine++) {
            cout << m_answer_horizontal[pixel(iLine, iCol)] << " ";
        }
        cout << "\n";
    }
}

void DD::printVerticalAnswers() {
    cout << "Printing vert :\n";
    for (int iLine=0; iLine<m_width; iLine++) {
        for (int iCol=0; iCol<m_height; iCol++) {
            cout << m_answer_vertical[pixel(iLine, iCol)] << " ";
        }
        cout << "\n";
    }
}


void DD::minSum(vector<vector<float> > &chain, vector<Label> &answer, int indexChain) {
    float minAct(Inf), marginal(Inf);
    Label labelOpt(0);
    int sizeChain = chain.size(), x(0), y(0);
    vector<vector<float> > forward(sizeChain,0), backward(sizeChain,0);
    for (int i=0; i<sizeChain ; i++) {
        forward[i].resize(m_nLabels);
        backward[i].resize(m_nLabels);
    }
    
    // Forward pass
    for (int node=1; node<sizeChain; node++) {
        for (Label labelNode=0; labelNode<m_nLabels; labelNode++) {
            minAct = Inf;
            for (Label  labelPreviousNode=0; labelPreviousNode<m_nLabels; labelPreviousNode++) {
                
                minAct = min(minAct, chain[node-1][labelPreviousNode] + forward[node-1][labelPreviousNode] + (float)m_V(labelNode, labelPreviousNode));
            }
            forward[node][labelNode] = minAct;
        }
    }
    
    // compute optimal labeling for the last node
    minAct = Inf; 
    for (int label=0; label<m_nLabels; label++) {
        marginal = chain[sizeChain-1][label] + forward[sizeChain-1][label];
        if (marginal < minAct) {
            minAct = marginal;
            labelOpt = label;
        }
    }
    if (indexChain < m_height) {
        y = indexChain;
        x = sizeChain - 1;
    }
    else {
        x = indexChain - m_height;
        y = sizeChain - 1;
    }
    answer[pixel(x, y)] = labelOpt;
    
    
    // backward pass + optimal labeling
    for (int node=sizeChain-2; node>=0; node--) {
        marginal = Inf;
        for (Label  labelNode=0; labelNode<m_nLabels; labelNode++) {
            minAct = Inf;
            for (Label  labelNextNode=0; labelNextNode<m_nLabels; labelNextNode++) {
                minAct = min(minAct, chain[node+1][labelNextNode] + backward[node+1][labelNextNode] + (float)m_V(labelNode, labelNextNode));
            }
            backward[node][labelNode] = minAct;
            minAct += forward[node][labelNode] + chain[node][labelNode];
            if (minAct<marginal) {
                marginal = minAct;
                labelOpt = labelNode;
            }
        }
        if (indexChain < m_height) {
            y = indexChain;
            x = node;
        }
        else {
            x = indexChain - m_height;
            y = node;
        }
        answer[pixel(x, y)] = labelOpt;
    }
    
} 


