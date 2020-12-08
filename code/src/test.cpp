#include <iostream>
#include <vector>
#include <map>
#include <string>
using namespace std;

// map<vector<int>, bool> taboList;

void print(const vector<int>& vecY, const string& str){
  cout << "vector " << str << " = {";
  for (auto& i : vecY)
    cout << i;
  cout << "}\n";
}

void taboSearch(const vector<int>& vecY){
  vector <int> vecYstar = vecY;
  int sz = vecY.size(), cnt = 0;
  print(vecY, "Y");
  for (int i = 0; i < sz; ++i)
    if (vecY[i] == 1)
      for (int j = i + 1; j < sz; ++j)
        if (vecY[j] == 1)
          for (int x = 0; x < sz; ++x)
            if (vecY[x] == 0)
              for (int y = x + 1; y < sz; ++y)
                if (vecY[y] == 0) {
                  // found
                  cnt++;
                  vecYstar = vecY;
                  vecYstar[i] = 0;
                  vecYstar[j] = 0;
                  vecYstar[x] = 1;
                  vecYstar[y] = 1;
                  print(vecYstar, "Y*(" + to_string(cnt) + ")");
                  // taboList[{i, j, x, y}] = true;
                }
}


int main(){
  vector<int> initialSol = {0, 0, 1, 1, 0, 0, 1, 0, 1, 1};
  taboSearch(initialSol);
}
