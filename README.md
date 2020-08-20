# 非線形シュレディンガー方程式を用いたビーム伝搬法（スカラ）をJuliaで。
スカラ波動方程式をビーム伝搬法で解く  
基本的にはYarivの論文をもとに解いていく。  
  
## 計算手法
3次元導波路のBeamPropagationMethod（BPM）を有限差分ADI法で解く。  
光波に対して遅延した形で屈折率が変調していくので、  
遅延時間ごとにセルを更新していく形にする。  
つまり、光波伝搬方向（Z方向）をZfまで計算し、  
Δt後にZf+ΔZまで計算するスキーム。  
Z方向は  
0→0+Δz  
0→Δz+Δz  
0→2Δz+Δz  
というように少しずつ伸ばしていく。  
  
ビーム伝搬法の解法には多様有り、
1. FFT-BPM（フーリエ変換BPM）  
    メリット：コーディングが簡単  
    デメリット：計算量、計算誤差が大きい
2. FE-BPM(有限要素法BPM)  
    メリット：複雑な形状にも対応  
    デメリット：計算量、コーディングがえぐい
3. FDTD－BPM(FDTD-BPM)  
    メリット：ほぼすべての状況に対応  
    デメリット：計算量がさらにえぐい。対称形などで使う。コーディングはもっとえぐい。
4. FD-BPM(有限差分法BPM)  
    メリット：コーディングが（やや）簡単  
    デメリット：矩形にのみ対応。

というようなことから、ここではFD-BPMを採用する。  

## 数式
