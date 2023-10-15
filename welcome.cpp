#include <iostream>
using namespace std;
void stamp(int* value, const char* var){
    cout << "inizializzo la variabile  " << var << " con " << *value << endl;
    cout << "inserisci un nuovo valore" ;
    cin >> *value;
    if (cin.fail ()) cout << "errore nella variabile in input" << endl;
    cout << "  la nuova variabile: " << *value << " con puntatore!" << endl;
}
int main(){
    int a=6;
    
    stamp(&a, "a");
    
    double vec[5]={0};
    for(int i=0; i<a; i++){
        cout << "iterazione n " << i << endl;
    }
    

    for (int i=0; i<5; i++){
        vec[i]=i*6;
     cout << "vec[" << i << "] = " << vec[i]<<endl;
    }
    return 0;
}
