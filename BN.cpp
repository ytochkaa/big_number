#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <algorithm>
#include <random>
#include <conio.h>
#include <vector>

using namespace std;


typedef unsigned int BASE;
#define BASE_SIZE (sizeof(BASE)*8) //16
typedef unsigned long long DBASE;

class BN {
    BASE* coef;
    int len;
    int maxlen;
public:
    explicit BN(int t = 1, int maxl = 1);
    BN(const BN&);
    ~BN() {
        if (coef) delete[]coef;
        coef = NULL;
    };
    BN& operator = (const BN&);
    void cout_16();
    void cin_16(const char*);
    void cin_10(int);

    int compare(const BN&) const;
    bool operator == (const BN&) const;
    bool operator != (const BN&) const;
    bool operator > (const BN&) const;
    bool operator >= (const BN&) const;
    bool operator < (const BN&) const;
    bool operator <= (const BN&) const;

    BN operator + (const BN&);
    BN operator += (const BN&);
    BN operator - (const BN&);
    BN operator -= (const BN&);
    BN operator * (const BASE);
    BN operator * (const BN&);
    BN operator *= (const BN&);
    BN operator + (const BASE);

    BN operator / (const BASE);
    BN operator % (const BASE);
    BN operator ^ (const BN&);
    void cout_10();
    void cin_10();

    BN operator /  (const BN&);
    BN operator %  (const BN&);
    void len_norm();
    friend void test();
    BN DICH(BASE);
    BN SQRT();
    BN SQRTB(int);
    BN pow_mod(const BN& exponent, const BN& mod) const;
    bool Miller_Rabin(int);
    void MPD();
    BN fast_square();
    BN Alway();
    void Ferma();
    void P0();
    void P1();
    int lg();

    
    void Gelfonds_algorithm(const BN&, int);
    void Group_Pol(BN& p, BN& g, BN& a, BN& x, BN& u, BN& v);
    void Ro_Pol(const BN&, const BN&);


    BN GCD(BN, BN);
    friend BN Fun(BN, BN);


};

BN::BN(int t, int maxl) {
    if (t == 1) {
        len = 1;
        maxlen = 1;
        coef = new BASE;
        coef[0] = 0;
    }
    else {
        if (t == 2) {
            maxlen = maxl;
            coef = new BASE[maxlen];
            len = 1;
            for (int i = 0; i < maxlen; i++) coef[i] = 0;
        }
        else {
            if (t == 3) {
                maxlen = maxl;
                coef = new BASE[maxlen];
                len = maxl;
                //maxl++;
                for (int i = 0; i < len; i++) {
                    coef[i] = rand();
                    if (sizeof(BASE) >= 4) {
                        coef[i] << 16;
                        coef[i] |= rand();
                    }
                }
                len_norm();
            }
            else {
                cout << "!!";
                exit(0);
            }
        }
    }
};

BN::BN(const BN& bn) {
    len = bn.len;
    maxlen = bn.maxlen;
    coef = new BASE[maxlen];
    for (int i = 0; i < maxlen; i++) {
        coef[i] = bn.coef[i];
    }
}

BN& BN::operator =(const BN& bn) {
    if (this != &bn) {
        delete[]coef;
        maxlen = bn.maxlen;
        len = bn.len;
        coef = new BASE[maxlen];
        for (int i = 0; i < maxlen; i++) {
            coef[i] = bn.coef[i];
        }
    }
    return *this;
};

void BN::cout_16() {
    int i = 0, tmp = 0;
    int k = BASE_SIZE - 4;
    char* s = new char[len * BASE_SIZE / 4];
    int j = 0;
    bool flag = true; int mask = 0;
    if (len == 1 && coef[0] == 0) cout << "0" << endl;
    else {
        for (j = len - 1; j >= 0;) {
            if (j == len - 1 && flag) {
                tmp = (coef[j] >> k) & (0xf);
                if (tmp >= 1 && tmp <= 9) {
                    s[i] = (char)(tmp + '0');
                    i++;
                    flag = false;
                }
                if (tmp >= 10 && tmp <= 15) {
                    s[i] = (char)(tmp - 10 + 'a');
                    i++;
                    flag = false;
                }
                k -= 4;
                if (k < 0) { k = BASE_SIZE - 4; j--; }
            }
            else {
                tmp = (coef[j] >> k) & (0xf);
                if (tmp >= 0 && tmp <= 9) {
                    s[i] = (char)(tmp + '0');
                }
                if (tmp >= 10 && tmp <= 15) {
                    s[i] = (char)(tmp - 10 + 'a');
                }
                i++; k -= 4;
                if (k < 0) { k = BASE_SIZE - 4; j--; }
            }
        }
        s[i] = '\0';
        puts(s);
    }
};

void BN::cin_16(const char* s) {
    int j = 0, k = 0, tmp = 0;
    BN TMP(2, ((strlen(s) - 1) / (BASE_SIZE / 4) + 1));
    TMP.len = TMP.maxlen;
    *this = TMP;
    for (int i = strlen(s) - 1; i >= 0; i--) {
        if (s[i] >= '0' && s[i] <= '9') {
            tmp = s[i] - '0';
        }
        else {
            if (s[i] >= 'a' && s[i] <= 'f') {
                tmp = s[i] - 'a' + 10;
            }
            else {
                if (s[i] >= 'A' && s[i] <= 'F') {
                    tmp = s[i] - 'A' + 10;
                }
                else { cout << "Error"; exit(1); }
            }
        }
        coef[j] |= tmp << (k * 4);
        k++;
        if (k >= BASE_SIZE / 4) { k = 0; j++; }
    }
    j = len - 1;
    while (coef[j] == 0 && j > 0) j--;
    len = j + 1;

};

int BN::compare(const BN& bn) const {
    if (len > bn.len) return 1;
    if (len < bn.len) return -1;
    for (int i = len - 1; i >= 0; i--) {
        if (coef[i] > bn.coef[i]) return 1;
        if (coef[i] < bn.coef[i]) return -1;
    }
    return 0;
};

bool BN::operator ==(const BN& bn) const {
    int tmp = this->compare(bn);
    if (tmp == 0) return 1; else return 0;
};

bool BN::operator !=(const BN& bn) const {
    int tmp = this->compare(bn);
    if (tmp == 1 || tmp == -1) return 1; else return 0;
};

bool BN::operator >(const BN& bn) const {
    int tmp = this->compare(bn);
    if (tmp == 1) return 1; else return 0;
};

bool BN::operator >=(const BN& bn) const {
    int tmp = this->compare(bn);
    if (tmp == 0 || tmp == 1) return 1; else return 0;
};

bool BN::operator <(const BN& bn) const {
    int tmp = this->compare(bn);
    if (tmp == -1) return 1; else return 0;
};

bool BN::operator <=(const BN& bn) const {
    int tmp = this->compare(bn);
    if (tmp == 0 || tmp == -1) return 1; else return 0;
};

BN BN:: operator + (const BN& bn) {
    int t, l;
    l = max(len, bn.len) + 1;
    t = min(len, bn.len);
    BN TMP(2, l);
    TMP.len = l;
    DBASE tmp;
    int j = 0, k = 0;
    while (j < t) {
        tmp = (DBASE)coef[j] + (DBASE)bn.coef[j] + (DBASE)k;
        TMP.coef[j] = (BASE)tmp;
        k = (BASE)(tmp >> BASE_SIZE);
        j++;
    }
    while (j < len) {
        tmp = (DBASE)coef[j] + (DBASE)k;
        TMP.coef[j] = (BASE)tmp;
        k = (BASE)(tmp >> BASE_SIZE);
        j++;
    }


    while (j < bn.len) {
        tmp = (DBASE)bn.coef[j] + (DBASE)k;
        TMP.coef[j] = (BASE)tmp;
        k = (BASE)(tmp >> BASE_SIZE);
        j++;
    }
    TMP.coef[j] = k;

    TMP.len_norm();
    return TMP;
};

BN BN:: operator += (const BN& bn) {
    *this = *this + bn;
    return *this;
};

BN BN:: operator - (const BN& bn) {
    if (bn.len > len) { cout << "Error.1 number must be bigger then 2"; exit(0); }
    DBASE tmp;
    BN TMP(2, len);
    TMP.len = len;
    int j = 0, k = 0;
    while (j < bn.len) {
        tmp = (DBASE)((DBASE)1 << BASE_SIZE) | (DBASE)coef[j]; 
        tmp = tmp - (DBASE)bn.coef[j] - (DBASE)k; 
        TMP.coef[j] = (BASE)tmp;
        k = !(tmp >> BASE_SIZE);
        j++;
    }
    while (j < len) {
        tmp = (DBASE)((DBASE)1 << BASE_SIZE) | (DBASE)coef[j];
        tmp -= (DBASE)k;
        TMP.coef[j] = (BASE)tmp;
        k = !(tmp >> BASE_SIZE);
        j++;
    }
    TMP.len_norm();
    return TMP;
};

BN BN:: operator -= (const BN& bn) {
    *this = *this - bn;
    return *this;
};

BN BN:: operator * (BASE v) {
    int j = 0;
    BASE k = 0;
    BN TMP(2, len + 1);
    TMP.len = TMP.maxlen;
    for (; j < len;) {
        DBASE tmp = (DBASE)coef[j] * (DBASE)v + (DBASE)k;
        TMP.coef[j] = (BASE)tmp;
        k = (BASE)(tmp >> BASE_SIZE);
        j++;
    }
    TMP.coef[j] = k;
    TMP.len = len + 1;
    TMP.len_norm();
    return TMP;
};

BN BN:: operator * (const BN& bn) {
    BN TMP(2, len + bn.len);
    TMP.len = TMP.maxlen;
    DBASE tmp;
    int j = 0;
    for (; j < bn.len;) {
        if (bn.coef[j] != 0) {
            int i = 0;
            BASE k = 0;
            for (; i < len;) {
                tmp = (DBASE)coef[i] * (DBASE)bn.coef[j] + (DBASE)TMP.coef[i + j] + (DBASE)k;
                TMP.coef[i + j] = (BASE)tmp;
                k = (BASE)(tmp >> BASE_SIZE);
                i++;
            }
            TMP.coef[j + len] = k;
        }
        j++;
    }
    TMP.len = len + bn.len;
    TMP.len_norm();
    return TMP;
};

BN BN:: operator *= (const BN& bn) {
    *this = *this * bn;
    return *this;
};

BN BN:: operator + (BASE v) {
    BN TMP(2, len + 1);
    TMP.len = len + 1;
    DBASE tmp;
    int j = 0, k = 0;
    tmp = (DBASE)coef[j] + (DBASE)v + (DBASE)k;
    TMP.coef[j] = (BASE)tmp;
    k = (BASE)(tmp >> BASE_SIZE);
    j++;

    while (j < len) {
        tmp = (DBASE)coef[j] + (DBASE)k;
        TMP.coef[j] = (BASE)tmp;
        k = (BASE)(tmp >> BASE_SIZE);
        j++;
    }
    TMP.coef[j] = k;
    TMP.len_norm();
    return TMP;
};

BN BN::operator / (BASE v) {
    if (v == 0) { cout << "Error! no / 0."; exit(0); }
    int j = 0; BASE r = 0;
    BN q(2, len);
    q.len = q.maxlen;
    DBASE tmp = 0;
    for (; j < len; j++) {
        tmp = ((DBASE)r << BASE_SIZE) + (DBASE)coef[len - 1 - j]; 
        q.coef[len - 1 - j] = (BASE)(tmp / (DBASE)v);
        r = (BASE)(tmp % (DBASE)v);
    }
    q.len = len;
    q.len_norm();
    return q;
};

BN BN:: operator % (BASE v) {
    if (v == 0) { cout << "Error! no / 0."; exit(0); }
    int j = 0; BASE r = 0;
    BN res(1);
    DBASE tmp = 0;
    for (; j < len; j++) {
        tmp = ((DBASE)r << BASE_SIZE) + (DBASE)coef[len - 1 - j];
        r = (BASE)(tmp % (DBASE)v);
    }
    res.coef[0] = r;
    return res;
};

void BN::cout_10() {
    string s;
    BN res = *this;
    BN zero(1);
    BN t(1, 1);

    while (res != zero) {
        t = res % 10;
        s.push_back(t.coef[0] + '0');
        res = res / 10;
        zero.len = res.len;
    }
    reverse(s.begin(), s.end());
    cout << s << " ";
}

void BN::cin_10() {
    char s[500];
    cin >> s;
    int j = 0;
    BASE tmp = 0;
    j = 0;
    BN TMP(2, (((strlen(s) - 1) / (BASE_SIZE / 4)) + 1));
    TMP.len = TMP.maxlen;
    if (s[0] == '-') {
        cout << "Error: Negative numbers are not allowed" << endl;
        // Ошибка: Отрицательные числа не разрешены
        exit(0);
    }
    for (; j < strlen(s); j++) {
        if (s[j] < '0' || s[j] > '9') { 
            cout << "Error: Invalid character in input" << endl; 
            // Ошибка: Недопустимый символ во входных данных
            return; 
        }
        tmp = (BASE)(s[j] - '0');
        TMP = TMP * 10 + tmp;
    }
    *this = TMP;
};

void BN::cin_10(int k) {
    string s = to_string(k);
    int j;
    BASE tmp = 0;
    BN TMP(2, (((s.length() - 1) / (BASE_SIZE / 4)) + 1));
    TMP.len = TMP.maxlen;
    if (s[0] == '-')
        exit(0);
    for (j = 0; j < s.length(); j++) {
        if (s[j] < '0' || s[j] > '9') {
            cout << "Invalid Input!";
            return;
        }
        tmp = (BASE)(s[j] - '0');
        TMP = TMP * 10 + tmp;
    }
    *this = TMP;
}

void BN::len_norm() {
    len = maxlen;
    while (len > 1 && coef[len - 1] == 0) {
        len--;
    }
}

BN BN:: operator / (const BN& bn) {

    if (bn.len == 1 && bn.coef[0] == 0) {
        throw invalid_argument("Error! No / 0.");
    }
    if (*this < bn) {
        BN q(1);
        return q;
    }

    if (bn.len == 1) {
        BN q;
        q = *this / bn.coef[0];
        return q;
    }

    int m = len - bn.len;
    DBASE b;
    b = (DBASE)((DBASE)1 << (DBASE)BASE_SIZE);
    DBASE d;
    d = (DBASE)(b / (DBASE)(bn.coef[bn.len - 1] + (BASE)1));
    int j = m;
    int k = 0;

    BN divisible;
    divisible = *this;
    divisible.maxlen++;
    divisible = divisible * d;
    BN divider;
    divider = bn;
    divider = divider * d;

    if (divisible.len == len) {
        divisible.coef[divisible.len++] = 0;
    }

    BN q(2, m + 1);
    q.len = m + 1;

    while (j >= 0) {
        DBASE q_ = ((((DBASE)divisible.coef[j + (DBASE)divider.len]) * (DBASE)(b)) + (DBASE)(divisible.coef[j + divider.len - 1])) / (DBASE)(divider.coef[divider.len - (DBASE)1]);
        DBASE r_ = ((((DBASE)divisible.coef[j + (DBASE)divider.len]) * (DBASE)(b)) + (DBASE)(divisible.coef[j + divider.len - 1])) % (DBASE)(divider.coef[divider.len - (DBASE)1]);

        if ((DBASE)q_ == (DBASE)b || (DBASE)(((DBASE)q_) * ((DBASE)divider.coef[divider.len - 2])) > (DBASE)((DBASE)b * ((DBASE)r_) + ((DBASE)divisible.coef[j + divider.len - 2]))) {
            q_--;
            r_ = (DBASE)r_ + ((DBASE)divider.coef[divider.len - 1]);
            if ((DBASE)r_ < (DBASE)b) {
                if ((DBASE)q_ == (DBASE)b || (DBASE)((DBASE)q_ * (DBASE)divider.coef[divider.len - 2]) > (DBASE)((DBASE)b * (DBASE)r_ + (DBASE)divisible.coef[j + divider.len - 2])) {
                    q_--;
                    r_ = (DBASE)r_ + (DBASE)divider.coef[divider.len - 1];
                }
            }
        }

        BN tmp(2, divider.len + 1);
        tmp.len = divider.len + 1;
        for (int i = 0; i < divider.len + 1; i++) {
            tmp.coef[i] = divisible.coef[j + i];
        }

        if (tmp < divider * (BASE)(q_)) {
            k = 1;
        }
        else k = 0;

        tmp = tmp - divider * (BASE)(q_);
        q.coef[j] = (BASE)(q_);

        if (k == 1) {
            q.coef[j]--;
            tmp = tmp + divider;
        }

        for (int i = 0; i < divider.len + 1; i++) {
            divisible.coef[j + i] = tmp.coef[i];
        }

        j--;
    }

    q.len_norm();

    return q;
}

BN BN:: operator % (const BN& bn) {

    if (bn.len == 1 && bn.coef[0] == 0) {
        throw invalid_argument("Error! No / 0.");
    }
    if (*this < bn) {
        return *this;
    }

    if (bn.len == 1) {
        BN r;
        r = *this % bn.coef[0];
        return r;
    }

    int m = len - bn.len;
    DBASE b;
    b = (DBASE)((DBASE)1 << (DBASE)BASE_SIZE);
    DBASE d;
    d = (DBASE)(b / (DBASE)(bn.coef[bn.len - 1] + (BASE)1));
    int j;
    j = m;
    int k = 0;

    BN divisible;
    divisible = *this;
    divisible.maxlen++;

    divisible = divisible * d;
    BN divider;
    divider = bn;
    divider = divider * d;

    if (divisible.len == len) {
        divisible.coef[divisible.len++] = 0;
    }

    while (j >= 0) {
        DBASE q_ = ((((DBASE)divisible.coef[j + (DBASE)divider.len]) * (DBASE)(b)) + (DBASE)(divisible.coef[j + divider.len - 1])) / (DBASE)(divider.coef[divider.len - (DBASE)1]);
        DBASE r_ = ((((DBASE)divisible.coef[j + (DBASE)divider.len]) * (DBASE)(b)) + (DBASE)(divisible.coef[j + divider.len - 1])) % (DBASE)(divider.coef[divider.len - (DBASE)1]);

        if ((DBASE)q_ == (DBASE)b || (DBASE)(((DBASE)q_) * ((DBASE)divider.coef[divider.len - 2])) > (DBASE)((DBASE)b * ((DBASE)r_) + ((DBASE)divisible.coef[j + divider.len - 2]))) {
            q_--;
            r_ = (DBASE)r_ + ((DBASE)divider.coef[divider.len - 1]);
            if ((DBASE)r_ < (DBASE)b) {
                if ((DBASE)q_ == (DBASE)b || (DBASE)((DBASE)q_ * (DBASE)divider.coef[divider.len - 2]) > (DBASE)((DBASE)b * (DBASE)r_ + (DBASE)divisible.coef[j + divider.len - 2])) {
                    q_--;
                    r_ = (DBASE)r_ + (DBASE)divider.coef[divider.len - 1];
                }
            }
        }

        BN tmp(2, divider.len + 1);
        tmp.len = divider.len + 1;
        for (int i = 0; i < divider.len + 1; i++) {
            tmp.coef[i] = divisible.coef[j + i];
        }

        if (tmp < divider * (BASE)(q_)) {
            k = 1;
        }
        else k = 0;

        tmp = tmp - divider * (BASE)(q_);

        if (k == 1) {

            tmp = tmp + divider;

        }

        for (int i = 0; i < divider.len + 1; i++) {
            divisible.coef[j + i] = tmp.coef[i];
        }

        j--;
    }
    divisible = divisible / d;

    divisible.len_norm();

    return divisible;
}

BN BN:: operator ^ (const BN& bn) {
    BN res;
    BN nul(1, 1);
    BN one; one.coef[0] = 1;
    BN two; two = one + one;
    res.coef[0] = 1;
    BN n = *this;
    BN tmp; tmp = bn;
    while (tmp > nul) {
        if (tmp % two == one) {
            res *= n;
        }
        n = n * n;
        tmp = tmp / two;
    }
    return res;
}

void test() {
    int M = 1000;
    int T = 1000;
    BN A, B, C, D;
    do
    {
        int n = rand() % M + 1;
        int m = rand() % M + 1;
        BN tmp(3, n);

        BN TMP(3, m);

        A = tmp;
        B = TMP;
        C = A / B;
        D = A % B;

        cout << "m = " << m << " ; ";
        cout << "n = " << n << " ; ";
        /*A.cout_10();
        B.cout_10();
        C.cout_10();
        D.cout_10();*/
        cout << "T = " << T << endl;

    } while (A == C * B + D && A - D == C * B && D < B && --T);

}

BN BN::DICH(BASE y) {
    BN z;
    z.coef[0] = 1;
    if (y == 0)
        return z;
    if (y == 1)
        return *this;

    BN q(1);
    q = *this;
    BASE mask = 1;
    BASE y1 = y;
    while (y1 != 0) {
        if (y & mask) z = z * q;
        mask = mask << 1;
        BN tmp; tmp = q;
        q = q * tmp;
        y1 = y1 >> 1;
    }
    return z;
}

BN BN::SQRT() {
    BN x;
    x = *this;
    BN x0;
    x0 = x;
    x = ((*this / x) + x) / 2;
    while (x < x0) {
        x0 = x;
        x = ((*this / x) + x) / 2;
    }
    return x0;
}

BN BN::fast_square() {
    BN x;
    x = *this;
    int n = x.len;
    BN y(2, 2 * n);
    y.len = y.maxlen;
    DBASE cu = 0;
    DBASE uv = 0;
    BASE v = 0;
    DBASE cuv = 0;
    for (int i = 0; i < n; i++) {
        uv = (DBASE)y.coef[2 * i] + (DBASE)x.coef[i] * (DBASE)x.coef[i];
        y.coef[2 * i] = (BASE)uv;
        cu = (uv >> BASE_SIZE);
        for (int j = i + 1; j < n; j++) {
            cuv = (DBASE)((DBASE)(BASE)y.coef[i + j] + (DBASE)((BASE)((DBASE)x.coef[i] * (DBASE)x.coef[j]) * 2) +
                (DBASE)(BASE)cu);
            v = (BASE)cuv;
            cu = (DBASE)((DBASE)((DBASE)(((DBASE)x.coef[i] * (DBASE)(x.coef[j])) >> BASE_SIZE) * (DBASE)2) +
                (DBASE)((DBASE)cu >> BASE_SIZE) + (DBASE)((DBASE)cuv >> BASE_SIZE));
            y.coef[i + j] = v;
        }
        y.coef[i + x.len] += (BASE)cu;
        y.coef[i + x.len + 1] += (BASE)(cu >> BASE_SIZE);
    }
    y.len = y.maxlen;
    for (int i = 2 * x.len - 1; i > -1; i--) {
        if (y.coef[i] == 0) {
            y.len--;
        }
        else {
            break;
        }
    }
    return y;
}

BN BN::pow_mod(const BN& exponent, const BN& mod) const {
    BN base;
    base = *this;
    BN result, zero, one, two;
    one.coef[0] = 1;
    result = one;
    two = one * 2;
    BN exp;
    exp = exponent;

    while (exp > zero) {
        if (exp % two == one) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp = exp / two;
    }

    return result;
}

BN BN::SQRTB(int n) {
    BN x, a;
    x = *this;
    BN x0;
    x0 = x;
    a = x;
    do {
        x0 = x;
        x = (x0 * (n - 1) + (a / x0.DICH(n - 1))) / n;
    } while (x0 > x);
    return x0;
}

bool BN::Miller_Rabin(int t) {
    BN n; n = *this;
    BN nul;
    BN one;
    one.coef[0] = 1;
    int pow_two = 0;
     
    BN a;
    BN two;
    two = one + one;
    a = two;

    BN n1;
    n1 = n - one;
    int s = 0;
    while ((n1 % two) == nul) {
        n1 = n1 / two;
        s++;
    }
    BN r;
    r = n1;
    BN y;
    BN tmp(3, len - 1);
    for (int i = 0; i < t; i++) {
        if (i == 1) a = two + one + two;
        if (i >= 2) {
            if (i % 2 == 0) a = a + (one * 2);
            else a = a + (one * 4);
        }
        if (a >= two || a <= n - two) {
            y = a.pow_mod(r, n);
            if (y != one and y != n - one) {
                int j = 1;
                while (j < s and y != n - one) {
                    y = (y * y) % n;
                    if (y == one) {
                        cout << "composite" << endl;
                        return false;
                    }
                    j++;
                }
                if (y != n - one) {
                    cout << "composite" << endl;
                    return 0;
                }
            }
        }
    }
    cout << "prime" << endl;
    return true;
}

void BN::MPD() {
    BN n(3, 3);

    BN nul;
    BN zero(1,0);
    BN one(1,1);
    BN two(1,2);
    one.coef[0] = 1;

    n = *this;
    
    if (n <= one) {
        cout << "n = 1 ";
        return ;
    }

    int pow_two = 0;
    while (n % 2 == nul) {
        n = n / 2;
        pow_two++;
    }

    if (n <= (one * 3)) {
        cout << '2' << "^" << pow_two << " * ";
        n.cout_10();
        return;
    }

    if (n.Miller_Rabin(10)) {
        n.cout_10();
        return;
    }

    BN d(1); // первый (пробный) делитель
    d = one * 3;
    int k = 1;
    bool flag = true;
    cout << endl;
    cout << "n = ";

    BN qrt;
    qrt = this->SQRT();
    BN q, r; // частное и остаток

    if (pow_two > 0 ){
        cout << '2' << "^" << pow_two << " * ";
    }

    while (n != one) {

        q = n / d;
        r = n % d;

        if (r == nul) {
            d.cout_10();
            n = q;
            if (q <= d) {
                n.cout_10();
                break;
            }
            cout << "* ";
            continue;
        }

        if (q > d) {
            k++;
            if (k == 2)
                d = d + (one * 2);
            else {
                if (k % 2 != 0)
                    d = d + (one * 2);
                else
                    d = d + (one * 4);
            }

            // больше условия
            if (d > qrt) {
                n.cout_10();
                return;
            }

        } else {
            n.cout_10();
            cout << endl;
            return;
        }
    }
    return;
}

BN BN::Alway() {
    //n = dk*qk + rk
    BN n(3, 3);
    BN nul, one, two, four;
    one.coef[0] = 1;
    two = one + one;
    four = one * 4;
    n = *this;

    while (n % 2 == nul) {
        n = n / 2;
    }

    if (n.Miller_Rabin(10)) {
        return n;
    }
    
    BN d = n.SQRTB(3) * two + one; //н
    BN s = n.SQRT(); //в

    cout << "[ ";
    d.cout_10();
    cout << "; ";
    s.cout_10();
    cout << "]";

    BN r1 = n % d;
    BN r2 = n % (d - two);
    BN q = ((n / (d - two)) - (n / d)) * 4;

    BN r;

    while (d <= s) {
        d += two;
        // r = 2*r1 - r2 + q < 0
        if ((r1 * 2 + q) < r2) {
            r = ((r1 * two) + d + q ) - r2;
            r2 = r1;
            r1 = r;
            q += four;
        }
        else {
            r = ((r1 * two) + q) - r2;
            r2 = r1;
            r1 = r;
        }

        while (r1 >= d) {
            r1 -= d;
            q -= four;
        }

        if (r1 == nul) {
            cout << "\nDiviser is found: ";
            d.cout_10();
            return d;
        }
    }
    cout << "\nNo dividers in this range!" << endl;
    return nul;
}
int BN::lg() {
    BASE r;
    r = coef[len - 1];
    BASE mask = 1;
    mask <<= 31;
    int d = 0;
    int i;
    for (i = 0; i < 32 && (r & mask) == 0; i++) {
        mask >>= 1;
    }
    d = 32 - (i + 1);
    d = d + 32 * (len - 1);
    return d;
}

BN Fun(BN x, BN n) {
    BN one, two;
    one.coef[0] = 1;
    two = one + one;
    x = (x.pow_mod(two, n) + one) % n;
    return x;
}

BN BN::GCD(BN a, BN b) {
    BN t, c, nul, one, two;
    one.coef[0] = 1;
    two = one + one;
    int k, i = 0, j = 0;
    while (a % two == nul) {
        a = a / two;
        i++;
    }
    while (b % two == nul) {
        b = b / two;
        j++;
    }
    if (i < j)
        k = i;
    else k = j;
    while (a != b) {
        while (a < b) {
            t = a;
            a = b;
            b = t;
        }
        c = a - b;
        while (c % two == nul) {
            c = c / two;
        }
        a = c;
    }
    a = two.DICH(k) * a;
    return a;
}


void BN::Ferma() {
    //n^2 = x^2 - y^2 = (x-y)*(x+y)

    BN n(3, 3);
    BN nul, one, two, x;
    one.coef[0] = 1;
    two = one + one;
    n = *this;

    // степень 2
    while (n % two == nul) {
        n = n / 2;
    }

    x = n.SQRTB(2);

    // Проверка, является ли n полным квадратом
    if (x.fast_square() == n) {
        cout << "a = b = ";
        x.cout_10();
        return;
    }

    BN z, y, a, b;

    unsigned long long count = 0;

    BN upper_bound = (n + one) / 2;

    // [ x ; (n+1)/2 ]
    cout << "Search range: [";
    x.cout_10();
    cout << "; ";
    upper_bound.cout_10();
    cout << "]" << endl;

    while (x != upper_bound) {
        x = x + one;
        count++;

        // урааааа мы снова достигли 10000 итераций
        if (count % 10000 == 0) {
            cout << "Iteration: " << count << ", current x: ";
            x.cout_10();
            cout << endl;
        }

        if (x == upper_bound) {
            cout << "n is prime! (Ferma)" << endl;
            cout << "Prime number: ";
            n.cout_10();
            cout << endl;
            return;

        } else {
            if (count == 1) {
                z = (x + n.SQRTB(2))*(x - n.SQRTB(2));
            } else {
                z = z + x * 2  - one;
            }
            //z = x.fast_square() - n;// вот тут вот мы улетаем непонятно куда))))
            y = z.SQRTB(2);

            if (y.fast_square() == z) {
                a = x + y;
                b = x - y;
                
                cout << "Found factors after " << count << " iterations:" << endl;
                cout << "a = ";
                a.cout_10();
                cout << endl;
                cout << "b = ";
                b.cout_10();
                cout << endl;
                cout << "Verification: a * b = ";
                (a * b).cout_10();
                cout << endl;

                return;
            }
        }
    }
    cout << "No factors found in the range." << endl;
    return;
}

//Метод основывается на свойствах отображений конечного множества в себя.
//Поскольку множество конечно,любая итеративная последовательность должна стать периодической
void BN::P0() {
    BN n;
    n = *this;
    BN a, b, d, t, nul, one, two;
    one.coef[0] = 1;
    two = one + one;

    d = one;

    // 1.
    a = two;
    b = two;
    
    while (d == one) {

        //2.
        a = Fun(a, n);
        b = Fun(Fun(b, n), n);

        //3. 
        if (a == b) {
            cout << "Without result";
            return;
        }
        
        //4.
        if (a < b){
            t = b - a;
        }else{
            t = a - b;
        }
        d = GCD(t, n);
    }
    
    cout << "Diviser is ";
    d.cout_10();
    cout << endl;
}
//исходя из некоторой границы гладкости B, найти делитель d числа n такой, что d – 1 – B-гладкое.
//B-гладкое, если все его простыеделители не превосходят B.
void BN::P1() {
    BN n;
    n = *this;
    BN a, q, d, one, tmp, B;
    int e;
    one.coef[0] = 1;

    a = one + one;

    int primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
                    109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229 };
    
    int primes_count = sizeof(primes) / sizeof(primes[0]);

//2.    
    for (int i = 0; i <= 35; i++) {
        d = GCD(a, n);
        if (d > one) {
            cout << "Diviser is ";
            d.cout_10();
            cout << endl;
            return;
        }
        
        BN a_current = a;
        bool success = false;
//3.        
        for (int j = 0; j < primes_count; j++) {
            q.cin_10(primes[j]);
            
            // e = [ln n / ln q ]
            int log_n = n.lg();
            int log_q = q.lg();

            if (log_q == 0) {
                log_q = 1; 
            }

            e = log_n / log_q;

            if (e < 1) {
                e = 1;
            }

            // Вычисляем q^e
            BN q_power = q.DICH(e);
            
            // a_current = a_current^(q^e) mod n
            a_current = a_current.pow_mod(q_power, n);
        }

//4:
        if (a_current == one) {
            a = a + one;
            continue;
        }
        tmp = a_current - one;
        d = GCD(tmp, n);
        
//5.        
        if (d > one && d < n) {
            cout << "Diviser is ";
            d.cout_10();
            cout << endl;
            return;  // Нашли нетривиальный делитель
        }
        a = a + one;
    }
    
    cout << "No divisor found" << endl;
}
/* gˣ ≡ a (mod p) */
/* x = log_g(a) */
void BN::Gelfonds_algorithm(const BN& P, int A) {
    BN one, p, g, n;
    one.coef[0] = 1;
    g = *this; p = P; n = p - one;
    
    // Конвертируем int A в BN a
    BN a;
    a.cin_10(A);
    
    // 1. h = [√n] + 1
    BN h = n.SQRT() + one;
    
    // 2. b = g^h mod p
    BN b = g.pow_mod(h, p);
    
    // Определяем размер h как целое число
    int h_int = 0;
    BN temp = h;
    BN zero; zero.coef[0] = 0;
    while (temp > zero) {
        h_int++;
        temp = temp - one;
    }
    
    // 3. Построение таблиц
    vector<BN> giant_steps; // b^u
    vector<BN> baby_steps;  // a * g^v
    
    // Шаги великана: b^u для u = 1, 2, ..., h
    BN u_val = one;
    for (int u = 1; u <= h_int; u++) {
        giant_steps.push_back(b.pow_mod(u_val, p));
        u_val = u_val + one;
    }
    
    // Шаги ребенка: a * g^v для v = 1, 2, ..., h  
    BN v_val = one;
    for (int v = 1; v <= h_int; v++) {
        baby_steps.push_back((a * g.pow_mod(v_val, p)) % p);
        v_val = v_val + one;
    }
    
    // 4. Поиск совпадения b^u = a * g^v
    for (int u = 0; u < h_int; u++) {
        for (int v = 0; v < h_int; v++) {
            if (giant_steps[u] == baby_steps[v]) {
                BN u_bn; u_bn.cin_10(u + 1);
                BN v_bn; v_bn.cin_10(v + 1);
                BN x = (h * u_bn - v_bn) % n;
                cout << "Result: ";
                x.cout_10();
                cout << endl;
                return;
            }
        }
    }
    
    cout << "No solution found" << endl;
}

void BN::Ro_Pol(const BN& P, const BN& A) {
    BN g, p, n, a, zero, one, two;
    zero.coef[0] = 0;
    one.coef[0] = 1;
    two = one + one;
    
    g = *this;
    p = P;
    n = p - one;
    a = A;
    
    // 1. Инициализация
    BN x1 = one, x2 = one;
    BN u1 = zero, u2 = zero, v1 = zero, v2 = zero;
    
    // 2-3. Поиск коллизии методом Флойда
    do {
        Group_Pol(p, g, a, x1, u1, v1);
        Group_Pol(p, g, a, x2, u2, v2);
        Group_Pol(p, g, a, x2, u2, v2);
    } while (x1 != x2);
    
    // 4. Вычисление r = v1 - v2 mod n
    BN r;
    if (v1 >= v2) {
        r = (v1 - v2) % n;
    } else {
        r = (n + v1 - v2) % n;
    }
    
    if (r == zero) {
        cout << "Error: v1 == v2, algorithm failed" << endl;
        return;
    }
    
    // Вычисление b = u2 - u1 mod n
    BN b;
    if (u2 >= u1) {
        b = (u2 - u1) % n;
    } else {
        b = (n + u2 - u1) % n;
    }
    
    // 5. Решение сравнения r*z ≡ b (mod n)
    BN d = GCD(r, n);
    
    if (d == one) {
        // Единственное решение
        // Находим z такое, что r*z ≡ b (mod n)
        BN z = zero;
        BN i = zero;
        while (i < n) {
            if ((r * i) % n == b) {
                z = i;
                break;
            }
            i = i + one;
        }
        
        cout << "Result: ";
        z.cout_10();
        cout << endl;
    } else {
        // Несколько решений
        BN r1 = r / d;
        BN b1 = b / d;
        BN n1 = n / d;
        
        // Находим z0 такое, что r1*z0 ≡ b1 (mod n1)
        BN z0 = zero;
        BN i = zero;
        while (i < n1) {
            if ((r1 * i) % n1 == b1) {
                z0 = i;
                break;
            }
            i = i + one;
        }
        
        // 6. Проверка всех решений z_i = z0 + k*n1, k=0..d-1
        BN k = zero;
        while (k < d) {
            BN z_candidate = z0 + k * n1;
            if (g.pow_mod(z_candidate, p) == a) {
                cout << "Result: ";
                z_candidate.cout_10();
                cout << endl;
                return;
            }
            k = k + one;
        }
        
        cout << "No valid solution found" << endl;
    }
}

void BN::Group_Pol(BN& p, BN& g, BN& a, BN& x, BN& u, BN& v) {
    BN one, two, three;
    one.coef[0] = 1;
    two = one + one;
    three = two + one;
    BN n = p - one; // порядок группы
    
    // Определяем, к какому подмножеству принадлежит x
    BN remainder = x % three;
    
    if (remainder == one) {
        // x ∈ S1: x = a * x
        x = (a * x) % p;
        v = (v + one) % n;
    }
    else if (remainder == two) {
        // x ∈ S2: x = x^2
        x = (x * x) % p;
        u = (u * two) % n;
        v = (v * two) % n;
    }
    else {
        // x ∈ S3: x = g * x
        x = (g * x) % p;
        u = (u + one) % n;
    }
}

int main() {
    BN g, p;
    int a_int;
    
    cout << "Format: g^x <=> a (mod p)\n" << endl;
    
    // Ввод основания g
    cout << "Enter base g (decimal): ";
    g.cin_10();
    
    // Ввод модуля p
    cout << "Enter modulus p (decimal): ";
    p.cin_10();
    
    // Ввод числа a как int
    cout << "Enter number a (decimal integer): ";
    cin >> a_int;
    
    cout << "\n ";
    g.cout_10();
    cout << "^ x <=> " << a_int << " ( mod ";
    p.cout_10();
    cout << ")" << endl;

    cout << "=== Gelfond's algorithm ===" << endl;
    g.Gelfonds_algorithm(p, a_int);
    
    // Для Ro_Pol конвертируем int в BN
    BN a_bn;
    a_bn.cin_10(a_int);
    
    cout << "=== Pollard's algorithm ===" << endl;
    g.Ro_Pol(p, a_bn);
    return 0;
}