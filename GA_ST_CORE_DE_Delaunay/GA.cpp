#include "GA.h"
#include "delaunator.hpp"
// Árvore Geradora Steady State
GA::GA(Populacao *pop, unsigned int qtdInd, float indCruz, float indMut, double fitAt) {
	this->populacao = pop;
	this->indCruzamento = indCruz;
	this->indMutacao = indMut;
	this->fitAtingir = fitAt;

	Individuo ind;
	for (unsigned int i = 0; i < qtdInd; i++) {
        this->populacao->addIndividuo(ind);   // Adiciona Indivíduo à população
    }
    unsigned int qtdPtFx = ind.getQtdPtFixo();
// Descobre triangulos de delaunay
    vector< double > coords;
    for (unsigned int i = 0; i < qtdPtFx; i++) {
        Vertice &iVert = ind.getVertice(i);
        coords.push_back(iVert.getX());
        coords.push_back(iVert.getY());
    }
    delaunator::Delaunator d(coords);

    // Triângulos de delaunay
    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        double tx0 = d.coords[2 * d.triangles[i]];
        double ty0 = d.coords[2 * d.triangles[i] + 1];
        double tx1 = d.coords[2 * d.triangles[i + 1]];
        double ty1 = d.coords[2 * d.triangles[i + 1] + 1];
        double tx2 = d.coords[2 * d.triangles[i + 2]];
        double ty2 = d.coords[2 * d.triangles[i + 2] + 1];

        double stX, stY;
        Aux::getPtTorricelli(tx0, ty0, tx1, ty1, tx2, ty2, stX, stY);

        Vertice ptSt(stX, stY);
        this->ptSt.push_back(ptSt);
    }
    unsigned int tamPtSt = this->ptSt.size();

// Constrói população
    unsigned int qtdMaxPtSt = ind.getQtdMaxSteiner();
    unsigned int dif = tamPtSt - qtdMaxPtSt;
    #pragma omp parallel for
    for (unsigned int i = 0; i < qtdInd; i++) {
        Individuo &ind = this->populacao->getIndividuo(i);
        vector<unsigned int> iPtSt;
        for (unsigned int j = 0; j < tamPtSt; j++ ) {
            iPtSt.push_back(j);
        }
        for (unsigned int j = 0; j < dif; j++ ) {
            unsigned int iRem = rand() % iPtSt.size();
            iPtSt.erase(iPtSt.begin() + iRem);
        }
        for (unsigned int j = 0; j < iPtSt.size(); j++ ) {
            ind.addPontoSteiner(this->ptSt.at(iPtSt.at(j)));
        }
        this->criarArvoreKruskal(ind);
	}
}

void GA::criarArvoreKruskal(Individuo &ind) {
// Gera MST com Kruskal, elimina invalidos e adiciona novos
    unsigned int qtdPtTotal = ind.getQtdPtTotal();   // Total
    unsigned int qtdPtFixo = ind.getQtdPtFixo();     // Inicio dos Steiners
    unsigned int qtdPtSt = ind.getQtdPtSteiner();    // Qtd de Pt Steiner

    vector< vector< double > > arestas;                         // Lista de arestas que podem ser inseridas (distância, vert1, vert2)
    for (unsigned int i = 0; i < qtdPtTotal; i++) {			    // Linha da meia matriz de adjacência
		for (unsigned int j = i + 1; j < qtdPtTotal; j++) {     // Coluna da meia matriz de adjancência
            vector<double> dist = {0, (double)i, (double)j};
            arestas.push_back(dist);
		}
    }
    unsigned int qtdArestas = arestas.size();
    #pragma omp parallel for
    for (unsigned int i = 0; i < qtdArestas; i++) {
        arestas[i][0] = ind.getDistancia(arestas[i][1], arestas[i][2]);
    }

    // Ordena lista de arestas em ordem crescente
    sort(arestas.begin(), arestas.end(), [](const vector< double >& a, const vector< double >& b){ return a[0] < b[0]; });

    // Cria arvore com base no algoritmo de Kruskal
    DisjointSets ds(qtdPtTotal);   // Cria conjuntos disjuntos

    // Percorre cromossomo
    for (unsigned int i = 0; i < qtdArestas; i++) {
        unsigned int u = arestas[i][1];
        unsigned int v = arestas[i][2];

        // Se um dos dois pontos for um Steiner com grau 3, não conecta
        if ((u >= qtdPtFixo && ind.getVertGrau(u) == 3) || (v >= qtdPtFixo && ind.getVertGrau(v) == 3)) {
            continue;
        }
        //if ((u < qtdPtFixo && ind.getVertGrau(u) == 1) || (v < qtdPtFixo && ind.getVertGrau(v) == 1)) {
            //continue;
        //}

        unsigned int set_u = ds.find(u);
        unsigned int set_v = ds.find(v);

        if (set_u != set_v) {                                   // Se não forma ciclo
            ds.merge(set_u, set_v);                             // Junta dois conjuntos
            ind.setGene(u, v, true);                            // Seta na matrix de adjacencia
        }
    }

	// Remove pontos Steiner com grau < 3
    unsigned int fim = qtdPtTotal;
	for (unsigned int i = 0; i < 3; i++) {      // Repete operação três vezes
        unsigned int iSt = qtdPtFixo;
        while (iSt < fim) {
            vector< unsigned int > adj;
            ind.getPtsAdjacentes(iSt, adj);
            if (adj.size() == 1) {
                ind.setGene(iSt, adj.at(0), false);       // Remove ligação do ponto Steiner
                ind.remPontoSteiner(iSt);
                fim--;
            } else if (adj.size() == 2) {
                ind.setGene(iSt, adj.at(0), false);       // Remove ligação do ponto Steiner
                ind.setGene(iSt, adj.at(1), false);       // Remove ligação do ponto Steiner
                ind.setGene(adj.at(0), adj.at(1), true);  // Liga dois pontos adjacentes ao Steiner
                ind.remPontoSteiner(iSt);
                fim--;
            } else {
                iSt++;
            }
        }
	}

	// Adiciona pontos Steiner e reposiciona
	this->criarArvoreSteiner(ind);
	this->reposicionaPontosSmith(ind);
}

void GA::iniciarGA(unsigned int numIteracoes, bool exibirInfo) {
    this->log = "";
    unsigned int qtdInd = this->populacao->getQtdIndividuos();
    double epsilon = 0.000005;

	unsigned int melhorInd = this->populacao->getDominante();             // Índice do melhor indivíduo da população
	Individuo melhor = this->populacao->getIndividuo(melhorInd);
    double melhorFit = melhor.getFitness();                               // Fitness do melhor indivíduo

    Individuo melhorGlob = melhor;
    double melhorFitGlob = melhorFit;

    //unsigned int qtdMaxStag = 4;//(numIteracoes * 2) / 100;        // 10% das iterações
    unsigned int qtdMaxStagSai = 10;//(numIteracoes * 20) / 100;    // 20% das iterações
    unsigned int qtdStag = 0;
    unsigned int qtdStagSai = 0;
    double melhorFitAnt = melhorFit;
    float f1 = this->indMutacao / 50;   // 0.5
    float f2 = this->indMutacao / 50;   // 0.5
    float cr = this->indCruzamento;     // 50%

    for (unsigned int iter = 1; iter <= numIteracoes; iter++) {
        if (qtdStag == qtdMaxStagSai) {
            break;
        }
        /*if (qtdStag == qtdMaxStag) {
            //for(unsigned int i = 0; i < qtdInd; i++) {
                //this->populacao->remIndividuo(0);
            //}
            //this->populacao->addIndividuo(melhor);
            this->populacao->remIndividuo(melhorInd);
            //for(unsigned int i = 1; i < qtdInd; i++) {
                Individuo individuo;
                this->criarArvorePrimRst(individuo);              // Cria árvore geradora mínima através do algoritmo de Prim para os pontos fixos
                this->criarArvoreSteiner(individuo);              // Adiciona pontos Steiner através de pontos fixos com grau > 2
                this->reposicionaPontosSmith(individuo);          // Ajusta posição
                this->populacao->addIndividuo(individuo);
            //}

            qtdStag = 0;
        }*/

        // Cria população temporária com indivíduos vazios
        Populacao popTemp;
        Individuo ind;
        unsigned int qtdMaxSt = ind.getQtdMaxSteiner();
        for (unsigned int i = 0; i < qtdMaxSt; i++) {
            Vertice ptSt = Vertice(0.0, 0.0);
            ind.addPontoSteiner(ptSt);
        }
        for (unsigned int iInd = 0; iInd < qtdInd; iInd++) {
            popTemp.addIndividuo(ind);
        }

        // Aplica cruzamento e mutação para criar os indivíduos da população temporária
        #pragma omp parallel for
        for (unsigned int iInd = 0; iInd < qtdInd; iInd++) {
            Individuo &filho = popTemp.getIndividuo(iInd);

            //setDE_rand_1(iInd, filho, f1, cr);                    // DE/rand/1
            //setDE_best_1(iInd, filho, f1, cr, melhorInd);         // DE/best/1
            setDE_rand_2(iInd, filho, f1, f2, cr);                // DE/rand/2
            //setDE_best_2(iInd, filho, f1, f2, cr, melhorInd);     // DE/best/2
            //setDE_curToRand_1(iInd, filho, f1, f2, cr);           // DE/ctor/1
            //setDE_curToBest_2(iInd, filho, f1, f2, cr, melhorInd); // DE/ctob/1

            this->criarArvoreKruskal(filho);                        // Cria topologia de árvore de Steiner
        }
        #pragma omp parallel for
        for (unsigned int iInd = 0; iInd < qtdInd; iInd++) {
            Individuo &paiCruz = this->populacao->getIndividuo(iInd);
            Individuo &filho = popTemp.getIndividuo(iInd);
            if (!isnan(filho.getFitness())) {
                if (filho.getFitness() < paiCruz.getFitness()) {
                    this->populacao->setIndividuo(iInd, filho);
                }
            }
        }

		melhorInd = this->populacao->getDominante();
		melhor = this->populacao->getIndividuo(melhorInd);
        melhorFit = melhor.getFitness();

        if (fabs(melhorFit - melhorFitAnt) < 0.0000000001) {
            qtdStag++;
            qtdStagSai++;
        } else {
            melhorFitAnt = melhorFit;
            qtdStag = 0;
            qtdStagSai = 0;
        }

        if (melhorFit < melhorFitGlob) {
            melhorFitGlob = melhorFit;
            melhorGlob = melhor;
        }

        // Armazena Log
        stringstream fitness;
        fitness << fixed << setprecision(9) << melhorFit;
        this->log += fitness.str() + "\n";
		//if (exibirInfo == true) {
            //double desvioPadrao = this->populacaoAGM->getDesvioPadrao();
            //cout << "Geração " << iter << "\tDesvio Padrão: " << desvioPadrao;
            //cout << "\tMelhor Fitness: " << this->populacaoAGM->getIndividuo(melhorInd).getFitness() << endl;
		//}

		//cout << setprecision(10) << melhorFit << " " << this->fitAtingir << endl;
		//if (fabs(melhorFit - this->fitAtingir) < epsilon) {
		if (melhorFit <= this->fitAtingir) {
            break;
		}
    }

    // Gera árvores de Steiner através das árvores geradoras
	this->reposicionaPontosSmith(melhorGlob, 0.000009);//, 0.000009
    this->populacao->addIndividuo(melhorGlob);
}

void GA::getSelecao(vector<unsigned int> &ind, unsigned int qtd, const unsigned int naoSel) {
    unsigned int qtdInd = this->populacao->getQtdIndividuos();

    vector<unsigned int> individuo;
    for (unsigned int i = 0; i < qtdInd; i++) {
        if (i !=  naoSel) {
            individuo.push_back(i);
        }
    }

    for (unsigned int i = 0; i < qtd; i++) {
        unsigned int rnd = rand() % individuo.size();
        ind.push_back(individuo.at(rnd));
        individuo.erase(individuo.begin() + rnd);
    }
}

void GA::setDE_rand_1(unsigned int iInd, Individuo &filho, float f1, float cr) {
    vector<unsigned int> iPaisMut;                                              // Indivíduos que participarão da mutação
    this->getSelecao(iPaisMut, 3, iInd);

    Individuo &paiCruz = this->populacao->getIndividuo(iInd);
    Individuo &paiMut1 = this->populacao->getIndividuo(iPaisMut.at(0));
    Individuo &paiMut2 = this->populacao->getIndividuo(iPaisMut.at(1));
    Individuo &paiMut3 = this->populacao->getIndividuo(iPaisMut.at(2));

    unsigned int qtdPtFx = filho.getQtdPtFixo();
    unsigned int qtdPtTot = filho.getQtdPtTotal();
    double xMin = filho.getXMin(); double xMax = filho.getXMax();
    double yMin = filho.getYMin(); double yMax = filho.getYMax();

    #pragma omp parallel for
    for (unsigned int iSt = qtdPtFx; iSt < qtdPtTot; iSt++) {                   // Percorre todos os pontos Steiner
        Vertice &i1 = paiMut1.getVertice(iSt);
        Vertice &i2 = paiMut2.getVertice(iSt);
        Vertice &i3 = paiMut3.getVertice(iSt);

        double x, y;
        if ((rand() % 100) < cr) {
            x = i1.getX() + (f1 * (i2.getX() - i3.getX()));
            //if (x < xMin) { x = xMin; } else if (x > xMax) { x = xMax; }
            if (x < xMin || x > xMax) { x = i1.getX(); }
            y = i1.getY() + (f1 * (i2.getY() - i3.getY()));
            //if (y < yMin) { y = yMin; } else if (y > yMax) { y = yMax; }
            if (y < yMin || y > yMax) { y = i1.getY(); }
        } else {
            x = paiCruz.getVertice(iSt).getX();
            y = paiCruz.getVertice(iSt).getY();
        }

        filho.getVertice(iSt).setX(x);
        filho.getVertice(iSt).setY(y);
    }
}

void GA::setDE_best_1(unsigned int iInd, Individuo &filho, float f1, float cr, unsigned int best) {
    vector<unsigned int> iPaisMut;                                              // Indivíduos que participarão da mutação
    this->getSelecao(iPaisMut, 2, iInd);

    Individuo &paiCruz = this->populacao->getIndividuo(iInd);
    Individuo &paiMut1 = this->populacao->getIndividuo(best);
    Individuo &paiMut2 = this->populacao->getIndividuo(iPaisMut.at(0));
    Individuo &paiMut3 = this->populacao->getIndividuo(iPaisMut.at(1));

    unsigned int qtdPtFx = filho.getQtdPtFixo();
    unsigned int qtdPtTot = filho.getQtdPtTotal();
    double xMin = filho.getXMin(); double xMax = filho.getXMax();
    double yMin = filho.getYMin(); double yMax = filho.getYMax();

    #pragma omp parallel for
    for (unsigned int iSt = qtdPtFx; iSt < qtdPtTot; iSt++) {                   // Percorre todos os pontos Steiner
        Vertice &i1 = paiMut1.getVertice(iSt);
        Vertice &i2 = paiMut2.getVertice(iSt);
        Vertice &i3 = paiMut3.getVertice(iSt);

        double x, y;
        if ((rand() % 100) < cr) {
            x = i1.getX() + (f1 * (i2.getX() - i3.getX()));
            //if (x < xMin) { x = xMin; } else if (x > xMax) { x = xMax; }
            if (x < xMin || x > xMax) { x = i1.getX(); }
            y = i1.getY() + (f1 * (i2.getY() - i3.getY()));
            //if (y < yMin) { y = yMin; } else if (y > yMax) { y = yMax; }
            if (y < yMin || y > yMax) { y = i1.getY(); }
        } else {
            x = paiCruz.getVertice(iSt).getX();
            y = paiCruz.getVertice(iSt).getY();
        }

        filho.getVertice(iSt).setX(x);
        filho.getVertice(iSt).setY(y);
    }
}

void GA::setDE_rand_2(unsigned int iInd, Individuo &filho, float f1, float f2, float cr) {
    vector<unsigned int> iPaisMut;                                              // Indivíduos que participarão da mutação
    this->getSelecao(iPaisMut, 5, iInd);

    Individuo &paiCruz = this->populacao->getIndividuo(iInd);
    Individuo &paiMut1 = this->populacao->getIndividuo(iPaisMut.at(0));
    Individuo &paiMut2 = this->populacao->getIndividuo(iPaisMut.at(1));
    Individuo &paiMut3 = this->populacao->getIndividuo(iPaisMut.at(2));
    Individuo &paiMut4 = this->populacao->getIndividuo(iPaisMut.at(3));
    Individuo &paiMut5 = this->populacao->getIndividuo(iPaisMut.at(4));

    unsigned int qtdPtFx = filho.getQtdPtFixo();
    unsigned int qtdPtTot = filho.getQtdPtTotal();
    double xMin = filho.getXMin(); double xMax = filho.getXMax();
    double yMin = filho.getYMin(); double yMax = filho.getYMax();

    #pragma omp parallel for
    for (unsigned int iSt = qtdPtFx; iSt < qtdPtTot; iSt++) {                   // Percorre todos os pontos Steiner
        Vertice &i1 = paiMut1.getVertice(iSt);
        Vertice &i2 = paiMut2.getVertice(iSt);
        Vertice &i3 = paiMut3.getVertice(iSt);
        Vertice &i4 = paiMut4.getVertice(iSt);
        Vertice &i5 = paiMut5.getVertice(iSt);

        double x, y;
        if ((rand() % 100) < cr) {
            x = i1.getX() + (f1 * (i2.getX() - i3.getX())) + (f2 * (i4.getX() - i5.getX()));
            //if (x < xMin) { x = xMin; } else if (x > xMax) { x = xMax; }
            if (x < xMin || x > xMax) { x = i1.getX(); }
            y = i1.getY() + (f1 * (i2.getY() - i3.getY())) + (f2 * (i4.getY() - i5.getY()));
            //if (y < yMin) { y = yMin; } else if (y > yMax) { y = yMax; }
            if (y < yMin || y > yMax) { y = i1.getY(); }
        } else {
            x = paiCruz.getVertice(iSt).getX();
            y = paiCruz.getVertice(iSt).getY();
        }

        filho.getVertice(iSt).setX(x);
        filho.getVertice(iSt).setY(y);
    }
}

void GA::setDE_best_2(unsigned int iInd, Individuo &filho, float f1, float f2, float cr, unsigned int best) {
    vector<unsigned int> iPaisMut;                                              // Indivíduos que participarão da mutação
    this->getSelecao(iPaisMut, 4, iInd);

    Individuo &paiCruz = this->populacao->getIndividuo(iInd);
    Individuo &paiMut1 = this->populacao->getIndividuo(best);
    Individuo &paiMut2 = this->populacao->getIndividuo(iPaisMut.at(0));
    Individuo &paiMut3 = this->populacao->getIndividuo(iPaisMut.at(1));
    Individuo &paiMut4 = this->populacao->getIndividuo(iPaisMut.at(2));
    Individuo &paiMut5 = this->populacao->getIndividuo(iPaisMut.at(3));

    unsigned int qtdPtFx = filho.getQtdPtFixo();
    unsigned int qtdPtTot = filho.getQtdPtTotal();
    double xMin = filho.getXMin(); double xMax = filho.getXMax();
    double yMin = filho.getYMin(); double yMax = filho.getYMax();

    #pragma omp parallel for
    for (unsigned int iSt = qtdPtFx; iSt < qtdPtTot; iSt++) {                   // Percorre todos os pontos Steiner
        Vertice &i1 = paiMut1.getVertice(iSt);
        Vertice &i2 = paiMut2.getVertice(iSt);
        Vertice &i3 = paiMut3.getVertice(iSt);
        Vertice &i4 = paiMut4.getVertice(iSt);
        Vertice &i5 = paiMut5.getVertice(iSt);

        double x, y;
        if ((rand() % 100) < cr) {
            x = i1.getX() + (f1 * (i2.getX() - i3.getX())) + (f2 * (i4.getX() - i5.getX()));
            //if (x < xMin) { x = xMin; } else if (x > xMax) { x = xMax; }
            if (x < xMin || x > xMax) { x = i1.getX(); }
            y = i1.getY() + (f1 * (i2.getY() - i3.getY())) + (f2 * (i4.getY() - i5.getY()));
            //if (y < yMin) { y = yMin; } else if (y > yMax) { y = yMax; }
            if (y < yMin || y > yMax) { y = i1.getY(); }
        } else {
            x = paiCruz.getVertice(iSt).getX();
            y = paiCruz.getVertice(iSt).getY();
        }

        filho.getVertice(iSt).setX(x);
        filho.getVertice(iSt).setY(y);
    }
}

void GA::setDE_curToRand_1(unsigned int iInd, Individuo &filho, float f1, float f2, float cr) {
    vector<unsigned int> iPaisMut;                                              // Indivíduos que participarão da mutação
    this->getSelecao(iPaisMut, 4, iInd);

    Individuo &paiCruz = this->populacao->getIndividuo(iInd);
    Individuo &paiMut1 = this->populacao->getIndividuo(iPaisMut.at(0));
    Individuo &paiMut2 = this->populacao->getIndividuo(iPaisMut.at(1));
    Individuo &paiMut3 = this->populacao->getIndividuo(iPaisMut.at(2));
    Individuo &paiMut4 = this->populacao->getIndividuo(iPaisMut.at(3));

    unsigned int qtdPtFx = filho.getQtdPtFixo();
    unsigned int qtdPtTot = filho.getQtdPtTotal();
    double xMin = filho.getXMin(); double xMax = filho.getXMax();
    double yMin = filho.getYMin(); double yMax = filho.getYMax();

    #pragma omp parallel for
    for (unsigned int iSt = qtdPtFx; iSt < qtdPtTot; iSt++) {                   // Percorre todos os pontos Steiner
        Vertice &i1 = paiMut1.getVertice(iSt);
        Vertice &i2 = paiMut2.getVertice(iSt);
        Vertice &i3 = paiMut3.getVertice(iSt);
        Vertice &i4 = paiMut4.getVertice(iSt);

        double x, y;
        if ((rand() % 100) < cr) {
            x = i1.getX() + (f1 * (i2.getX() - i1.getX())) + (f2 * (i3.getX() - i4.getX()));
            //if (x < xMin) { x = xMin; } else if (x > xMax) { x = xMax; }
            if (x < xMin || x > xMax) { x = i1.getX(); }
            y = i1.getY() + (f1 * (i2.getY() - i1.getY())) + (f2 * (i3.getY() - i4.getY()));
            //if (y < yMin) { y = yMin; } else if (y > yMax) { y = yMax; }
            if (y < yMin || y > yMax) { y = i1.getY(); }
        } else {
            x = paiCruz.getVertice(iSt).getX();
            y = paiCruz.getVertice(iSt).getY();
        }

        filho.getVertice(iSt).setX(x);
        filho.getVertice(iSt).setY(y);
    }
}

void GA::setDE_curToBest_2(unsigned int iInd, Individuo &filho, float f1, float f2, float cr, unsigned int best) {
    vector<unsigned int> iPaisMut;                                              // Indivíduos que participarão da mutação
    this->getSelecao(iPaisMut, 3, iInd);

    Individuo &paiCruz = this->populacao->getIndividuo(iInd);
    Individuo &paiMut1 = this->populacao->getIndividuo(best);
    Individuo &paiMut2 = this->populacao->getIndividuo(iPaisMut.at(0));
    Individuo &paiMut3 = this->populacao->getIndividuo(iPaisMut.at(1));
    Individuo &paiMut4 = this->populacao->getIndividuo(iPaisMut.at(2));

    unsigned int qtdPtFx = filho.getQtdPtFixo();
    unsigned int qtdPtTot = filho.getQtdPtTotal();
    double xMin = filho.getXMin(); double xMax = filho.getXMax();
    double yMin = filho.getYMin(); double yMax = filho.getYMax();

    #pragma omp parallel for
    for (unsigned int iSt = qtdPtFx; iSt < qtdPtTot; iSt++) {                   // Percorre todos os pontos Steiner
        Vertice &i1 = paiMut1.getVertice(iSt);
        Vertice &i2 = paiMut2.getVertice(iSt);
        Vertice &i3 = paiMut3.getVertice(iSt);
        Vertice &i4 = paiMut4.getVertice(iSt);

        double x, y;
        if ((rand() % 100) < cr) {
            x = i2.getX() + (f1 * (i1.getX() - i2.getX())) + (f2 * (i3.getX() - i4.getX()));
            //if (x < xMin) { x = xMin; } else if (x > xMax) { x = xMax; }
            if (x < xMin || x > xMax) { x = i1.getX(); }
            y = i2.getY() + (f1 * (i1.getY() - i2.getY())) + (f2 * (i3.getY() - i4.getY()));
            //if (y < yMin) { y = yMin; } else if (y > yMax) { y = yMax; }
            if (y < yMin || y > yMax) { y = i1.getY(); }
        } else {
            x = paiCruz.getVertice(iSt).getX();
            y = paiCruz.getVertice(iSt).getY();
        }

        filho.getVertice(iSt).setX(x);
        filho.getVertice(iSt).setY(y);
    }
}

unsigned int GA::addPontoSteiner(Individuo &ind, unsigned int v1, unsigned int v2, unsigned int v3) {		// Private
	// Insere um ponto Steiner no ponto de Torricelli do triangulo formado por 3 vértices conectados
    // Triângulo é formado por: v1 conectado a v2 e a v3
    // Retorna índice do vértice inserido

    double steinerX, steinerY;       // Posição dos pontos Steiner
    double x1 = ind.getVertice(v1).getX(); double y1 = ind.getVertice(v1).getY();
    double x2 = ind.getVertice(v2).getX(); double y2 = ind.getVertice(v2).getY();
    double x3 = ind.getVertice(v3).getX(); double y3 = ind.getVertice(v3).getY();

    // Cria ponto Steiner
    Aux::getPtTorricelli(x1, y1, x2, y2, x3, y3, steinerX, steinerY);
    Vertice pontoSteiner(steinerX, steinerY);
    ind.addPontoSteiner(pontoSteiner);
    unsigned int ultimoPt = ind.getQtdPtTotal() - 1;			// Indice do ponto Steiner inserido

    // Retira a arestas que liga (v1 a v2) e (v1 a v3), e os liga ao ponto steiner
    ind.setGene(v1, v2, false);
    ind.setGene(v1, v3, false);
    ind.setGene(v1, ultimoPt, true);
    ind.setGene(v2, ultimoPt, true);
    ind.setGene(v3, ultimoPt, true);

    return ultimoPt;    // Retorna índice
}

void GA::remPontoSteiner(Individuo &ind, unsigned int s1, unsigned int v1, unsigned int v2, unsigned int v3) {		// Private
	// Remove arestas que ligam ponto Steiner e liga v1-v2, v1-v3
	// Remove ponto Steiner
    ind.setGene(s1, v1, false);
    ind.setGene(s1, v2, false);
    ind.setGene(s1, v3, false);
    ind.setGene(v1, v2, true);
    ind.setGene(v1, v3, true);
    ind.remPontoSteiner(s1);
}

void GA::criarArvoreSteiner(Individuo &ind) {
    unsigned int qtdPtFixo = ind.getQtdPtFixo();

    // Insere a quantidade máxima de vértices Steiner criando uma "Árvore de Steiner cheia"
    vector<unsigned int> adjacente;
	for (unsigned int i = 0; i < qtdPtFixo; i++) {
		adjacente.clear();
        ind.getPtsAdjacentes(i, adjacente);
        unsigned int tamAdj = adjacente.size();

        while (tamAdj > 1) {        // Se possuir grau > 1
            // Se possuir grau 2
            unsigned int iVert = 0;
            unsigned int jVert = 1;
            if (tamAdj > 2) {   // Se o grau for maior que 2
                // Busca pelo conjunto de vértices adjacentes com menor ângulo
                double menorAngulo = 6.28319;                // 360° em radianos
                for (unsigned int iAdj = 0; iAdj < tamAdj - 1; iAdj++) {
                    for (unsigned int jAdj = iAdj + 1; jAdj < tamAdj; jAdj++) {
                        double angulo = ind.getAngulo(i, adjacente[iAdj], adjacente[jAdj]);
                        if (angulo < menorAngulo) {
                            iVert = iAdj;
                            jVert = jAdj;
                            menorAngulo = angulo;
                        }
                    }
                }
            }
            // Adiciona ponto Steiner
            unsigned int ultimoPt = this->addPontoSteiner(ind, i, adjacente[iVert], adjacente[jVert]);
            adjacente.erase(adjacente.begin() + iVert); jVert--;
            adjacente.erase(adjacente.begin() + jVert);
            adjacente.push_back(ultimoPt);
            tamAdj--;
        }
	}
}

void GA::reposicionaPontosSmith(Individuo &ind, double prec) {
    unsigned int qtdPtFixo = ind.getQtdPtFixo();            // Indice inicial dos pontos Steiner
    unsigned int qtdPtSteiner = ind.getQtdPtSteiner();      // Quantidade de pontos Steiner

    // Cria lista de adjacência e distância dos pontos Steiner
    vector< vector< unsigned int > > lstAdj(qtdPtSteiner, vector< unsigned int >(3));   // Lista de Adjacência dos pontos Steiner
    vector< vector< double > > lstDist(qtdPtSteiner, vector< double >(3));              // Lista de Distâncias dos pontos Steiner

    #pragma omp parallel for
    for (unsigned int i = 0; i < qtdPtSteiner; i++) {
        unsigned int indicePt = qtdPtFixo + i;
        vector<unsigned int> adj;
        ind.getPtsAdjacentes(indicePt, adj);
        if (adj.size() == 3) {
            lstAdj[i] = adj;
        } else {
            lstAdj[i].erase(lstAdj[i].begin());
        }
    }

    double q = this->length(ind, lstAdj, lstDist);
    double r = this->error(ind, lstAdj, lstDist);
    do {
        this->optimize(ind, lstAdj, lstDist, prec * (r / qtdPtFixo));
        q = this->length(ind, lstAdj, lstDist);
        r = this->error(ind, lstAdj, lstDist);
    } while ( r > (q * prec));                          // Usar 0.01 para performance e 0.00001 para precisão

    ind.setFitness(q);
}

void GA::optimize(Individuo &ind, vector< vector<unsigned int> > &lstAdj, vector< vector< double > > &lstDist, double tol) {
    unsigned int qtdPtFixo = ind.getQtdPtFixo();            // Indice inicial dos pontos Steiner
    unsigned int qtdPtSteiner = ind.getQtdPtSteiner();      // Quantidade de pontos Steiner

    // 1°: Computa B array, C array, e valences. Configura o leafQ.
    double B[qtdPtSteiner][3];                      // Adjacentes
    double C[qtdPtSteiner][2];                      // X e Y
    unsigned int val[qtdPtSteiner], leafQ[qtdPtSteiner], lqp;
	unsigned int eqnstack[qtdPtSteiner], eqp;

	#pragma omp parallel for
    for (unsigned int i = 0; i < qtdPtSteiner; i++) {
        if (lstAdj[i].size() != 3) { continue; }    // Se o ponto Steiner não possuir grau 3
        double q[3];                                // Pontos Adjacentes
		q[0] = 1.0 / (lstDist[i][0] + tol);
		q[1] = 1.0 / (lstDist[i][1] + tol);
		q[2] = 1.0 / (lstDist[i][2] + tol);
		double soma = q[0] + q[1] + q[2];
		q[0] /= soma;
		q[1] /= soma;
		q[2] /= soma;

		val[i] = 0;
		B[i][0] = B[i][1] = B[i][2] = 0.0;          // Pontos Adjacentes
		C[i][0] = C[i][1] = 0.0;                    // X e Y

		for (unsigned int j = 0; j < 3; j++) {      // Adjacentes
            unsigned int adj = lstAdj[i][j];
            if (adj >= qtdPtFixo) { 	            // Se adj for um ponto Steiner
                val[i]++;
                B[i][j] = q[j];
            } else {                                // Se for um ponto fixo
                Vertice &ptAdj = ind.getVertice(adj);
                C[i][0] += ptAdj.getX() * q[j];
                C[i][1] += ptAdj.getY() * q[j];
            }
		}
    }
    lqp = 0;
    for (unsigned int i = 0; i < qtdPtSteiner; i++) {
        if (lstAdj[i].size() != 3) { continue; }    // Se o ponto Steiner não possuir grau 3
        if (val[i] <= 1) {                          // Se for adjacente a dois ou mais Steiner
			leafQ[lqp] = i;                         // Insere folha em leafQ
			lqp++;                                  // Incrementa contador de folhas
		}
    }
    // Configurar equações para resolvê-las.
    // 2°: eliminar folhas
    eqp = 0;
    while (lqp > 1) {
		lqp--;
		unsigned int i = leafQ[lqp];                // Indice do vértice para listas
		val[i]--;
		unsigned int i2 = i + qtdPtFixo;            // Indice do ponto Steiner

		// Elimina folha i
		eqnstack[eqp] = i;                          // Insere i na pilha
		eqp++;

		unsigned int j;
		for (j = 0; j < 3; j++) {
            if (B[i][j] != 0.0 && lstAdj[i][j] >= qtdPtFixo) break;     // Adjacente é j
		}
        if (j == 3) {
            for (j = 0; j < 3; j++) {
                if (lstAdj[i][j] >= qtdPtFixo) break;                   // Adjacente é j
            }
        }

		double q0 = B[i][j];
		unsigned int j2 = lstAdj[i][j] - qtdPtFixo;                     // Adjacente é j

		val[j2]--;
		if (val[j2] == 1) {
			leafQ[lqp] = j2;
			lqp++;
		}

		// Nova folha?
		unsigned int m;
		for (m = 0; m < 3; m++) {
            if (lstAdj[j2][m] == i2) break;
		}

		double q1 = B[j2][m];
		B[j2][m] = 0.0;
		double t = 1.0 - (q1 * q0);
		t = 1.0 / t;
		for (m = 0; m < 3; m++) {
            B[j2][m] *= t;
		}
		for (m = 0; m < 2; m++) {
			C[j2][m] += q1 * C[i][m];
			C[j2][m] *= t;
		}
	}

    // 3°: Resolve 1-vertex tree
	unsigned int i = leafQ[0];
	unsigned int i2 = i + qtdPtFixo;
    ind.getVertice(i2).setX(C[i][0]);
    ind.getVertice(i2).setY(C[i][1]);

    // 4°: Resolve novamente
	while (eqp > 0) {
		eqp--;
		unsigned int i = eqnstack[eqp];
		unsigned int i2 = i + qtdPtFixo;
		unsigned int j;
		for (j = 0; j < 3; j++) {
            if (B[i][j] != 0.0 && lstAdj[i][j] >= qtdPtFixo) break; // Adjacente é j
		}
		if (j == 3) {
            for (j = 0; j < 3; j++) {
                if (lstAdj[i][j] >= qtdPtFixo) break;               // Adjacente é j
            }
        }
		double q0 = B[i][j];
		j = lstAdj[i][j];                           // Adjacente é j
        Vertice &vertJ = ind.getVertice(j);
        double x = C[i][0] + q0 * vertJ.getX();
        double y = C[i][1] + q0 * vertJ.getY();
        ind.getVertice(i2).setX(x);
        ind.getVertice(i2).setY(y);
	}
}

double GA::error(Individuo &ind, vector< vector<unsigned int> > &lstAdj, vector< vector< double > > &lstDist) {
    unsigned int qtdPtFixo = ind.getQtdPtFixo();
    unsigned int qtdPtSteiner = ind.getQtdPtSteiner();

    double efig = 0.0;                                      // Erro
    #pragma omp parallel for reduction(+: efig)
    for (unsigned int i = 0; i < qtdPtSteiner; i++) {
        if (lstAdj[i].size() != 3) { continue; }            // Se o ponto Steiner não possuir grau 3
        unsigned int indicePtSt = qtdPtFixo + i;            // Indice do ponto Steiner
        Vertice &ptSt = ind.getVertice(indicePtSt);
        Vertice &adj0 = ind.getVertice(lstAdj[i][0]);
        Vertice &adj1 = ind.getVertice(lstAdj[i][1]);
        Vertice &adj2 = ind.getVertice(lstAdj[i][2]);

        double d12, d01, d02, r, s, t;
        d12 = d01 = d02 = 0.0;
        for (unsigned int eixo = 0; eixo < 2; eixo++) {     // X e Y
            t = ptSt.getPos(eixo);
            r = adj0.getPos(eixo) - t;
            s = adj1.getPos(eixo) - t;
            t = adj2.getPos(eixo) - t;
            d12 += s * t;
            d01 += r * s;
            d02 += r * t;
        }
        // Somente angulos < 120 causam erro
		t = d12 + d12 + (lstDist[i][1] * lstDist[i][2]);
		if (t > 0.0) efig += t;
		t = d01 + d01 + (lstDist[i][0] * lstDist[i][1]);
		if (t > 0.0) efig += t;
		t = d02 + d02 + (lstDist[i][0] * lstDist[i][2]);
		if (t > 0.0) efig += t;
    }
    return sqrt(efig);
}

double GA::length(Individuo &ind, vector< vector<unsigned int> > &lstAdj, vector< vector< double > > &lstDist) {
    unsigned int qtdPtFixo = ind.getQtdPtFixo();
    unsigned int qtdPtSteiner = ind.getQtdPtSteiner();

    double distTot = 0.0;                                           // Distância total
    #pragma omp parallel for reduction(+: distTot)
    for (unsigned int i = 0; i < qtdPtSteiner; i++) {
        if (lstAdj[i].size() != 3) { continue; }                    // Se o ponto Steiner não possuir grau 3
        unsigned int indicePtSt = qtdPtFixo + i;                    // indice do Ponto Steiner
        for (unsigned int j = 0; j < 3; j++) {
            unsigned int indicePtAdj = lstAdj[i][j];                // Adjacente ao ponto Steiner

            if (indicePtAdj < indicePtSt) {
                double dist = ind.getDistancia(indicePtSt, indicePtAdj);
                distTot += dist;
                lstDist[i][j] = dist;                               // Insere dist em lista de distância

                if (indicePtAdj >= qtdPtFixo) {                     // Se adj for um ponto Steiner
                    indicePtAdj -= qtdPtFixo;                       // Indice no vetor de adjacencia e distância
                    for (unsigned int k = 0; k < 3; k++) {          // Busca ponto Steiner em adj
                        if (lstAdj[indicePtAdj][k] == indicePtSt) {
                            lstDist[indicePtAdj][k] = dist;         // Insere dist em lista de distância
                            break;
                        }
                    }
                }
            }
        }
    }
    return distTot;
}

string GA::getLog() {
    return this->log;
}
