
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cmath>

using namespace std;

/// General helpers
float **newMat(int rows, int cols);
template <class T>
	T *newVect(size_t dim);
void displayMat(size_t dim, float **A);
void displayVec(size_t dim, float *v);

class SolveCholeski{
public:
	SolveCholeski(string data);    //data est nom du fichier contenant les données de la matrice
	void displayResult();
    void decomposeMat();            //retourne une matrice triangulaire sup B
    void transposeMat(float **mat); //transposer une matrice
	void solveTriangSup();       //resolution de matrice triangulaire superieure
    void solveTriangInf();       //resolution de matrice triangulaire inferieure
    
	size_t  getdim(){ return dim;}
	float** getMat(){ return A;}
    float** getB(){return B;}
    float** getBt(){return Bt;}
	float*  getRhs(){ return b;}
private:
	size_t dim;
	float  **A;		    // matrice du probleme A.x=b
    float **B;          // matrice triangulaire sup d udecomposition
	float **Bt;         // matrice transposé de B
    float  *b, *x, *y;	// second membre et inconnu du probleme
};

SolveCholeski::SolveCholeski(string data){
/// ouverture du fichier de donnees
	size_t i(0), j(0);
    ifstream fichier(data, ios::in);
    if(fichier){
        fichier >> dim;
/// allocation de la matrice, du second membre, B, Bt, y et de la solution
        A = newMat(dim, dim);
        b = newVect<float> (dim);   //vecteur formé du second membre
        B = newMat(dim,dim);        //matrice de decomposition, triang sup
        Bt = newMat(dim,dim);       //matrice transposee de B
        x = newVect<float> (dim);   //solution du systeme
        y = newVect<float> (dim);   //vecteur tel que B.y = b

/// remplissage des donnees
        for(i=0; i<dim; i++){
            for(j=0; j<dim; j++){
                fichier >> A[i][j]; 
            }
        }
        for(i=0; i<dim; i++){
            fichier >> b[i];
		}
        for(i=0; i<dim; i++){       //initialisation des valeurs des vecteurs
            x[i] = 0;
            y[i] = 0;
		}
        fichier.close();            //On ferme le fichier
    }
    else{
        cout << "On ne trouve pas des données..." << endl;  //message d'erreur
    }
}

void SolveCholeski::displayResult(){
    cout << "\nVecteur y:" << endl;
    for(size_t i=0; i<(size_t)dim; i++)
        cout << "y" << i+1 << " = " << y[i] << endl;
    cout << "\nLa solution" << endl;
    for(size_t i=0; i<dim; i++)
        cout << "x" << i+1 << " = " << x[i] << endl;
}


void SolveCholeski::decomposeMat(){        //trouver B,Bt tq A=B.Bt et A.x=b
    for(int i=0; i<(int)dim;i++){
        for(int j=0; j<=i;j++){
            float sum = 0;                //somme des B[i][k].B[j][k]
            if(i!=j){
                for(int k=0;k<=j-1;k++){
                    sum+=(B[i][k]*B[j][k]);
                }
                B[i][j] = (1/B[j][j]*(A[i][j]-sum));
            }
            else if (i==j){       //i==j
                for(int k=0;k<=i-1;k++){
                    sum += pow(B[i][k],2);
                }
                B[i][i]= sqrt(A[i][i]- sum);
            }
        }
    }
    transposeMat(B);                             //remplir Bt avec le transposé de B
}

void SolveCholeski::transposeMat(float **mat){   //remplir Bt avec le transposé de la matrice mat
    for(int i=0;i<(int)dim;i++){
        for(int j=0;j<(int)dim;j++){
            Bt[i][j] = mat[j][i];
        }
    }
}

void SolveCholeski::solveTriangSup(){
    float s(0);
    int i(0), j(0);
    for(i=dim-1; i>=0; i--){    
        for(j=i+1, s=0; j<int(dim); j++)
            s += (Bt[i][j]*x[j]);
        x[i] = (y[i]-s)/Bt[i][i];
    }
}
void SolveCholeski::solveTriangInf(){
    float s(0);
    int i(0), j(0);
    for(i=0; i<(int)dim; i++){  
        for(j=0, s=0; j<int(dim); j++)
            s += (B[i][j]*y[j]);
        y[i] = (b[i]-s)/B[i][i];
    }
}

int main(){

    cout << "Resolution d'un systeme lineaire A.x=b par la méthode de Cholesky" << endl;
    
/// Les données dans data.txt 

//data.txt form:
/*
4
1 2 4 7
2 13 23 38
4 23 77 122
7 38 122 294
2
10
36
-43
*/
    SolveCholeski solver("data.txt");
    cout << "La matrice A : " << endl;
    displayMat(solver.getdim(), solver.getMat());
    cout << "Le second membre b : " << endl;
    displayVec(solver.getdim(), solver.getRhs());
    
/// Calculs
	solver.decomposeMat();     
    //display
    cout << "\n Decomposistion de la matrice A tel que A = B.Bt" << endl;
    cout << "\nLa matrice B :" << endl;
    displayMat(solver.getdim(),solver.getB());

    cout << "\nLa matrice Bt :" << endl;
    displayMat(solver.getdim(),solver.getBt());
    solver.solveTriangInf();
    solver.solveTriangSup();
/// Resultats
	solver.displayResult();

    return 0;
}

void displayMat(size_t dim, float **A){
	for(size_t i=0; i<dim; i++){
		if(i==0)cout << "[";
		else    cout << "[";
		for(size_t j=0; j<dim; j++){
			if(j<dim-1)cout << setw(6) << setprecision(5) << A[i][j] << "  ";
			else cout << setw(8) << setprecision(5) << A[i][j];
		}
		if(i<dim-1) cout << "]" << endl;
		else 		cout << "]" << endl;
	}
}

void displayVec(size_t dim, float *v){
	for(size_t i=0; i<dim; i++){
		cout << " [" << v[i] << "]" << endl;
	}
}

template <class T>
T *newVect(size_t dim){
    T *v(NULL);
    v = new T[dim];
    if(v==NULL){
        cout << "Erreur lors de la création du vecteur" << endl;
    }
    return v;
}

float **newMat(int rows, int cols){
    float **mat(NULL);
    mat = newVect<float *>(rows);

    for(int i=0; i<rows; i++)
        mat[i] = newVect<float>(cols);
    return mat;
}