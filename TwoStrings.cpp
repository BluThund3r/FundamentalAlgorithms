#include <bits/stdc++.h>
using namespace std;

class TwoStrings{
    char str1[256], str2[256];
public:
    TwoStrings(char s1[], char s2[]) {
        strcpy(str1, s1);
        strcpy(str2, s2);
    }

    void setStr1(char s[]) { strcpy(str1, s); }
    void setStr2(char s[]) { strcpy(str2, s); }
    int getEditDistance();
    
    //SUBSIR -> nu este contiguu
    int getLCSubsequenceLength();

    // SUBSECVENTA -> este contiguu
    int getLCSubstringLength();
};

int TwoStrings::getEditDistance() {
    int n = strlen(str1), m = strlen(str2);
    vector<vector<int>> d(n + 1, vector<int>(m + 1));

    for(int i = 0; i <= n; ++ i)
        d[i][0] = i;

    for(int i = 0; i <= m; ++ i)
        d[0][i] = i;
    
    for(int i = 1; i <= n; ++ i)
        for(int j = 1; j <= m; ++ j)
            d[i][j] = min(d[i - 1][j] + 1, min(d[i][j - 1] + 1, d[i - 1][j - 1] + (str1[i - 1] != str2[j - 1])));
        
    return d[n][m];
}

// SUBSECVENTA -> este contiguu
int TwoStrings::getLCSubstringLength() {
    int n = strlen(str1), m = strlen(str2);
    vector<vector<int>> d(n + 1, vector<int>(m + 1, 0));

    int maxInMatrix = 0;
    for(int i = 1; i <= n; ++ i)
        for(int j = 1; j <= m; ++ j){
            if(str1[i - 1] == str2[j - 1])
                d[i][j] = d[i - 1][j - 1] + 1;
            maxInMatrix = max(maxInMatrix, d[i][j]);
        }
             
    return maxInMatrix;
}

//SUBSIR -> nu este contiguu
int TwoStrings::getLCSubsequenceLength() {
    int n = strlen(str1), m = strlen(str2);
    vector<vector<int>> d(n + 1, vector<int>(m + 1, 0)), Max(n + 1, vector<int>(m + 1, 0));

    for(int i = 1; i <= n; ++ i)
        for(int j = 1; j <= m; ++ j) {
            if(str1[i - 1] == str2[j - 1])
                d[i][j] = 1 + Max[i - 1][j - 1];
            Max[i][j] = max(Max[i - 1][j], max(Max[i][j - 1], d[i][j]));
        }

    return d[n][m];
}

int main() {

    return 0;
}