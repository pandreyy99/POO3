//
// Created by Andrei on 15-May-18.
//

#ifndef POO3_POLINOM_H
#define POO3_POLINOM_H

#include <iostream>
#include <cmath>

using namespace std;

template<class T>
class polinom {
private:
    int grad; /// Gradul polinomului
    T *p; /// Vector de coeficienti care va fi alocat dinamic

public:
    /** Default constructor */
    polinom();

    /** Constructor cu initializare valori(o constanta)*/
    polinom(const int, int);

    /** Constructor cu parametri */
    polinom(int grd, T *coef);

    /** Default destructor */
    virtual ~polinom();

    /** Copy constructor
     *  \param other Object to copy from
     */
    polinom(const polinom<T> &other);

    //int& polinom_set( int , int * ) ;

    /** Calcularea valorii intr-un punct x */
    T pdx(T x);

    /** Aflarea gradului polinomului */
    int get();

    /** Adaugarea unui element(x) de rang i */
    void ad(int i, T x);

    /**Eliminarea elementului de rang i */
    void el(const int i);

    /** Elementul de rang i */
    T &operator[](int);

    /** Supraincarcarea operatorului = */
    polinom<T> &operator=(polinom<T>);

    /** Supraincarcarea operatorului + */
    template<class U>
    friend polinom<U> operator+(const polinom<U> &pol1, const polinom<U> &pol2);

    /** Supraincarcarea operatorului * (intre 2 polinoame) */
    template<class U>
    friend polinom<U> operator*(const polinom<U> &pol1, const polinom<U> &pol2);

    /** Supraincarcarea operatorului * (cu scalar) */
    /** Cand scalarul e al doilea factor */
    template<class U>
    friend polinom<U> operator*(const polinom<U> &pol, const U &a);

    /** Cand scalarul e primul factor */
    template<class U>
    friend polinom<U> operator*(const U &a, const polinom<U> &pol);

    /** Supraincarcarea operatorului / (intre 2 polinoame) */
    template<class U>
    friend polinom<U> operator/(const polinom<U> &pol1, const polinom<U> &pol2);

    /** Supraincarcarea operatorului de citire */
    template<class U>
    friend istream &operator>>(istream &x, polinom<U> &pol);

    /** Supraincarcarea operatorului de afisare */
    template<class U>
    friend ostream &operator<<(ostream &x, const polinom<U> &pol);
};

/**
 * definire metode
 */

template<class T>
polinom<T>::polinom() {
    ///ctor fara parametri
    grad = 0;
    p = new T[1];
    p[0] = (T) 0;
}

template<class T>
polinom<T>::polinom(const int c, int grd) {
    ///ctor cu coeficienti egali cu o constanta si cu un grad dat(primit ca parametru)
    grad = grd;
    p = new T[grad + 1];
    for (int i = 0; i <= grad; i++)
        p[i] = (T) c;
}

template<class T>
polinom<T>::polinom(int grd, T *coef) {
    ///ctor cu parametri(un grad si un vector de coef)
    grad = grd;
    p = new T[grad + 1];
    for (int i = 0; i <= grad; i++)
        p[i] = coef[i];
}

template<class T>
polinom<T>::~polinom() {
    ///dtor
    grad = 0;
    delete[]p;
}

template<class T>
polinom<T>::polinom(const polinom<T> &other) {
    ///copy ctor
    grad = other.grad;
    p = new T[grad + 1];
    for (int i = 0; i <= grad; i++)
        p[i] = other.p[i];
    /**if(this != &other){
            delete []p ;
        grad = other.grad ;
        p = new int[grad + 1] ;
        for( int i = 0 ; i <= grad ; i++ )
            p[i] = other.p[i] ;
    }*/
}

template<class T>
T polinom<T>::pdx(T x) {
    ///valoarea polinomului in punctul x
    int i;
    T s = 0;
    for (i = 0; i <= grad; i++)
        s += p[i] * pow(x, i);
    return s;
}

template<class T>
int polinom<T>::get() {
    ///aflarea gradului unui polinom , precum si actualizarea polinomului si a gradului in cazul in care nu mai sunt de actualitate
    if (p[grad] == (T) 0) {
        T *temp;
        int i, ok = 0;
        for (i = grad - 1; i >= 0 && ok == 0; i--)
            if (p[i] != (T) 0) {
                grad = i;
                ok = 1;
            }
        temp = new T[grad + 1];
        for (i = 0; i <= grad; i++)
            temp[i] = p[i];
        delete[]p;
        p = new T[grad + 1];
        for (i = 0; i <= grad; i++)
            p[i] = temp[i];
        delete[]temp;
    }
    return grad;
}

template<class T>
void polinom<T>::ad(int i, T x) {
    /// adaugarea valorii x pe pozitia i
    /// daca pozitiei i ii corespunde deja un coeficient , nu se va face nimic
    if (i <= grad) return;
    /// in caz contrar , ne folosim de un temporar
    T *temp;
    int j;
    /// realocam vectorul de coef dupa gradul actual
    temp = new T[i + 1];
    /// punem valoarea x pe pozitia i
    temp[i] = x;
    /// setam valoarea 0 coeficientilor de pe poz grad_vechi+1...i-1
    for (j = grad + 1; j < i; j++)
        temp[j] = (T) 0;
    /// copiem valorile vechi ale polinomului
    for (j = 0; j <= grad; j++)
        temp[j] = this->p[j];
    /// dezalocam zona veche de memorie
    delete[]this->p;
    this->p = new T[i + 1];
    this->grad = i;
    /// copiem noile valori ale polinomului
    for (j = 0; j <= i; j++)
        p[j] = temp[j];
    delete[]temp;
}

template<class T>
void polinom<T>::el(const int i) {
    /// daca pozitia e mai mare ca si gradul , nu facem nimic
    if (i > grad) return;
    /// daca exista ,si e diferita de poz coef dominant , setam coef aferent ei la 0
    if (i < grad) p[i] = 0;
    /// daca vrem sa eliminam coef dominant si implicit sa scadem gradul polinomului
    if (i == grad) {
        /// ne folosim de un temp
        T *temp;
        int j;
        /// alocam spatiul necesar noului polinom
        temp = new T[i];
        /// copiem valorile vechiului polinom pana la penultima poz
        for (j = 0; j < i; j++)
            temp[j] = p[j];
        /// dezalocam zona de memorie corespunzatoare coef vechiului polinom
        delete[]this->p;
        /// scadem gradul
        this->grad = i - 1;
        this->p = new T[i];
        /// copiem noile valori
        for (j = 0; j <= i; j++)
            p[j] = temp[j];
        delete[]temp;
    }
}

template<class T>
T &polinom<T>::operator[](int i) {
    T *a = new T;
    *a = (T) 0;
    /// daca exista elementul p[i] , il returnam
    if (i >= 0 && i <= grad)
        return p[i];
        /// altfel , returnam valoarea 0
    else return *a;
}

template<class U>
istream &operator>>(istream &x, polinom<U> &pol) {
    /// supraincarcarea op >> pt citirea unui polimom
    delete[]pol.p;\
    U local;
    x >> local;
    pol.grad = local;
    cout << pol.grad << endl;
    pol.p = new U[pol.grad + 1];
    for (int i = pol.grad; i >= 0; i--)
        x >> pol.p[i];
    /// returnam fluxul
    return x;
}

template<class U>
ostream &operator<<(ostream &x, const polinom<U> &pol) {
    /// supraincarcarea op << pt afisarea unui polinom in forma corespunzatoare
    if ((pol.grad == 0) && (pol.p[pol.grad] == 0)) x << "0 ";
    else {
        if (pol.p[pol.grad] != 1)
            x << pol.p[pol.grad] << " * X^" << pol.grad;
        else x << "X^" << pol.grad;
        for (int i = (pol.grad) - 1; i >= 0; i--) {
            if (i != 0 && i != 1) {
                if (pol.p[i] < 0) {
                    if (pol.p[i] != -1)
                        x << " - " << (-1) * pol.p[i] << " * X^" << i;
                    else x << " - " << "X^" << i;

                }
                if (pol.p[i] > 0) {
                    if (pol.p[i] != 1)
                        x << " + " << pol.p[i] << " * X^" << i;
                    else x << " + " << "X^" << i;
                }
            }
            if (i == 1) {
                if (pol.p[1] < 0) x << " - " << (-1) * pol.p[1] << " * X";
                if (pol.p[1] > 0) {
                    if (pol.p[1] != 1) x << " + " << pol.p[1] << " * X";
                    else x << " + " << " X";
                }
            }
            if (i == 0) {
                if (pol.p[0] < 0) x << " - " << (-1) * pol.p[0];
                if (pol.p[0] > 0) x << " + " << pol.p[0];
            }
        }
    }
    return x;
}

template<class T>
polinom<T> &polinom<T>::operator=(polinom<T> pol) {
    T *temp = p;
    p = new T[pol.grad + 1];
    grad = pol.grad;
    for (int i = 0; i <= grad; i++)
        p[i] = pol.p[i];
    delete[]temp;
    return *this;
}

template<class U>
polinom<U> operator+(const polinom<U> &pol1, const polinom<U> &pol2) {
    polinom<U> pol;
    int i;
    /// gradul polinomului rezultat va fi = gradul maxim dintre cele 2 polinoame
    pol.grad = max(pol1.grad, pol2.grad);
    pol.p = new U[pol.grad + 1];
    /// daca primul grad va fi >= cel de-al doilea grad , adunam coef de pe poz 0,gradul comun si doar ii atribuim pe cei din primul polinom pentru celelalte
    if (pol1.grad >= pol2.grad) {
        for (i = 0; i <= pol2.grad; i++)
            pol.p[i] = pol1.p[i] + pol2.p[i];
        for (i = pol2.grad + 1; i <= pol1.grad; i++)
            pol.p[i] = pol1.p[i];
    }
        /// analog daca cel de-al doilea grad e mai mare decat primul
    else {
        for (i = 0; i <= pol1.grad; i++)
            pol.p[i] = pol1.p[i] + pol2.p[i];
        for (i = pol1.grad + 1; i <= pol2.grad; i++)
            pol.p[i] = pol2.p[i];
    }
    return pol;
}

template<class U>
polinom<U> operator*(const polinom<U> &pol1, const polinom<U> &pol2) {
    polinom<U> pol(0, pol1.grad + pol2.grad);
    int i, j;
    /// gradul polinomului rezultat va fi suma gradelor polinoamelor factor
    pol.grad = pol1.grad + pol2.grad;
    for (i = 0; i <= pol1.grad; i++)
        for (j = 0; j <= pol2.grad; j++)
            pol.p[i + j] += pol1.p[i] * pol2.p[j];
    return pol;
}

template<class U>
polinom<U> operator*(const polinom<U> &pol, const U &a) {
    /// inmultirea unui polinom cu un scalar(scalarul fiind al 2 - lea factor)
    polinom<U> pol_temp(0, pol.grad);
    int j;
    U i, s = 0;
    for (j = 0; j <= pol.grad; j++) {
        i = pol.p[j];
        pol_temp.p[j] += pol.p[j] * a + s;
        s += i;
    }
    return pol_temp;
}

template<class U>
polinom<U> operator*(const U &a, const polinom<U> &pol) {
    /// cand scalarul e primul factor
    polinom<U> pol_temp(0, pol.grad);
    int i;
    for (i = 0; i <= pol.grad; i++)
        pol_temp.p[i] = a * pol.p[i];
    return pol_temp;
}

template<class U>
polinom<U> operator/(const polinom<U> &pol1, const polinom<U> &pol2) {
    /// impartirea a 2 polinoame
    /// gradul polinomului rezultat va fi egal cu abs(diferenta celor 2 grade)
    U a = -1;
    int grd = abs(pol1.grad - pol2.grad);
    polinom<U> pol(0, grd), temp(0, grd), temp1(pol1), temp2(pol2);
    /// daca gradul impartitorului e mai mare decat cel al deimpartitului returnam polinomul nul
    if (pol1.grad < pol2.grad) {
        pol.grad = 0;
        pol.p[0] = 0;
        return pol;
    }
    ///else
    pol.grad = temp1.grad - temp2.grad;
    temp.grad = temp1.grad - temp2.grad;
    /// atata timp cat gradul deimpartitului e >= decat gradul impartitorului(impartirea e posibila)
    ///aplicam algoritmul lui euclid
    while (temp1.get() >= temp2.get()) {
        temp.p[temp1.grad - temp2.grad] = temp1.p[temp1.grad] / temp2.p[temp2.grad];
        temp.grad = temp1.grad - temp2.grad;
        pol = pol + temp;
        /// inmultim cu -1 polinomul temp(in caz contrar era necesara supraincarcarea op -)
        temp = temp * a;
        /// folosim teorema impartirii cu rest
        temp = temp * temp2;
        /// actualizam polinomul temp
        temp.grad = temp.get();
        /// daca coef catului partial e 0 , ne oprim
        if (temp.p[temp.grad] == 0) temp1.grad = 0;
        else {
            /// actualizam deimpartitul
            temp1 = temp1 + temp;
            /// setam valoarea coef ce a dat rangul catului temporar la 0(eroare cauzata de convertirea la int a catului)
            temp1.p[temp1.grad] = 0;
        }
    }
    return pol;
}

#endif //POO3_POLINOM_H