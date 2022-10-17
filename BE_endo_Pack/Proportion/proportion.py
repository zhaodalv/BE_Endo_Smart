import pickle
import numpy as np
import math
import pickle
import time
import itertools
import os
import copy
import pandas as pd


def solve(a, b, c, d): 
    t1 = b * (d - 1) / a #a
    t2 = d + b / a - c / a * (d - 1) ##b
    t3 = -c / a #c
    y = (-t2 + math.sqrt(t2 * t2 - 4 * t1 * t3)) / (2 * t1)
    x = (c - b * y) / a  ##
    return x, y

class BayesianRes():
    def __init__(self, prob):
        self.prob = prob
        self.length = prob.shape[0]

    def infer(self, obs):
        pre = 0
        ret = 1.0
        obs = np.array(obs, dtype=np.int64)
        for i in range(len(obs)):
            if obs[i] == 1:
                #print("1:",self.prob[i][pre])
                ret *= self.prob[i][pre]
            else:
                #print('0:',1-self.prob[i][pre])
                ret *= 1-self.prob[i][pre]
            pre = obs[i]
            #print(pre)
        #print(ret)
        return ret

    def inferall(self):
        return_dict = {}
        for i in itertools.product("01", repeat=self.length):
            return_dict[i] = self.infer(list(i))

        return return_dict


class BayesianNetwork():
    def __init__(self, score, positions=None, metric="cross") -> None:
        if positions is not None:
            self.positions = positions
        else:
            self.positions = list(range(20))
        #self.pos, self.neg = countfrequency(data, positions)

        self.score = score
        # self.solve = solve

        if metric == "corr":
            self.score[np.isnan(self.score)] = 0
            self.solve = solvecorr
        elif metric == "cross":
            self.score[np.isnan(self.score)] = 1
            self.solve = solve
        else:
            raise NotImplementedError

    def fit(self, positions, values):
        subset = []
        for i in range(len(positions)):
            if positions[i]:
                if i in self.positions:
                    subset.append(i)
        subset.sort()
        pro = np.zeros((len(subset), 2))

        if len(subset) == 0:
            pro = np.zeros((1, 2))
            pro[0] = 1.0
            return BayesianRes(pro)

        pro[0][1] = values[subset[0]]
        pro[0][0] = pro[0][1]

        

        for i in range(1, len(subset)):
            x, y = self.solve(
                values[subset[i - 1]], 1 - values[subset[i - 1]],
                values[subset[i]], self.score[self.positions.index(
                    subset[i - 1])][self.positions.index(subset[i])])
            #print(values[subset[i-1]]*x+(1-values[subset[i-1]])*y, values[subset[i]])
            # print((values[subset[i-1]]*x-values[subset[i-1]]*values[subset[i]])/math.sqrt(values[subset[i-1]]*(1-values[subset[i-1]])*values[subset[i]]*(1-values[subset[i]])),
            #  self.score[self.positions.index(subset[i-1])][self.positions.index(subset[i])])
            pro[i][1] = x
            pro[i][0] = y
        res = BayesianRes(pro)
        return res


class ProportionResults():
    def __init__(self, seq, predict, predictpos, solver,editor):
        self.mapping = {}
        self.predictpos = predictpos
        self.editor=editor
        for i in range(len(predictpos)):
            if predictpos[i]:
                self.mapping[len(self.mapping)] = i
        self.seq = seq
        self.res = solver.fit(predictpos, predict) ##predict 20bp label
        #print(self.mapping)
        #print(self.predictpos)
        self.infer = None
        
    def get_seq(self,t):
        change="T"
        if self.editor=='CBE':
            change='T'
        else:
            change='G'
        sgRNA=self.seq[10:30] ##40bp 10~30 is sgRNA
        sgRNA_list=list(sgRNA)
        out_seq=[]
        for x,y in zip(sgRNA_list,t):
            if y:
                out_seq.append(change)
            else:
                out_seq.append(x)
        return ''.join(out_seq)
    
    def inferall(self):
        if self.infer is None:
            self.infer = {}
        else:
            return self.infer
        ret = self.res.inferall()
        #print(ret)
        for i, v in ret.items():
            t = copy.deepcopy(self.predictpos)
            
            for j in range(len(self.mapping)):
                #print(i,j)
                if i[j] == "0":
                    t[self.mapping[j]] = False
            #print(t)
            out_seq=self.get_seq(t)
            self.infer[out_seq]=v
            #self.infer[str(t)] = v
        return self.infer

class Proportion_prediction:
    def __init__(self,editor,pre_bayes):
        self.Base_bayes = BayesianNetwork(pre_bayes)
        self.ET=editor
        
    def prediction(self,Seq40s,sgRNA,label,eff):
        Raw_predict_proportion=[]
        for ID,sgRNA,label,eff in zip(Seq40s,sgRNA,label,eff):
            result=ProportionResults(ID,list(eff),list(label),self.Base_bayes,self.ET)
            infer_proportion=result.inferall()
            infer_proportion_df=pd.DataFrame.from_dict(infer_proportion,orient='index')
            infer_proportion_df.columns=['Predicted_eff']
            infer_proportion_df.sort_values('Predicted_eff',ascending=False,inplace=True)
            infer_proportion_df['Seq_40']=ID
            infer_proportion_df['target']=sgRNA
            infer_proportion_df['outcomes']=infer_proportion_df.index
            Raw_predict_proportion.append(infer_proportion_df)
            del(result)
        Raw_predict_proportion_df=pd.concat(Raw_predict_proportion)
        return Raw_predict_proportion_df[['Seq_40','target','outcomes','Predicted_eff']]
