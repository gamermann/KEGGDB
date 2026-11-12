# Copyright (C) 2026 Daniel Gamermann
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from time import time
from functools import reduce
from .paths import *
import os
import re
import pkgutil


#####################################
###### Internal stuff
####################################


def _replacenn(stri: str, nn: int=13, mm: int=3) -> str:
    """
    
    Replaces generic indication of a number.

    :param stri: str: 
    :param nn: int:  (Default value = 13)
    :param mm: int:  (Default value = 3)

    """
    if type(stri)!=str:
        return stri
    stri = re.sub(r'(\d)([a-z])', '\\1*\\2', stri)
    strin = []
    for ele in stri:
        if not ele.isalpha():
            strin.append(ele)
        elif ele=='n':
            strin.append(str(nn))
        else:
            strin.append(str(mm))
    return eval("".join(strin))

def _compChemNames(name1: str, name2:str, penalty: float=10, thres: float=.05) -> float:
    """
    
    Compares the names of compounds.
    
    :param name1: str: 
    :param name2: str: 
    :param penalty: float:  (Default value = 10)
    :param thres: float:  (Default value = .05)

    """
    name1 = name1.lower()
    name2 = name2.lower()
    l1 = len(name1)
    l2 = len(name2)
    diff = abs(l2-l1)/(l2+l1)
    if l1==l2 or diff<thres:
        poin = 0.
        ll = min(l2, l1)
        return sum([name1[ii]!=name2[ii] for ii in range(ll)])/ll
    elif l1<l2:
        if name1 in name2:
            return (l2-l1)/l2
        else:
            return penalty
    else:
        if name2 in name1:
            return (l1-l2)/l1
        else:
            return penalty



def _colapseMultipleDics(dics: list[chemDic[str, float]], nums: list[int]|None=None) -> chemDic[str, float]:
    """
    Transforms many scoring dictionaries into one.
    
    :param dics: A list with the scoring dictionaries to be combined.
    :type dics: list[chemDic[str, float]]
    :param nums: A list of numbers by which each dictionary (score within) should be multiplied. Default is None (a list of ones).
    :type nums: list[int]|None

    :returns: The combined dictionary. 
    :rtype: chemDic[str, float]
    """
    if nums is None:
        nums = [1 for ele in dics]
    if not dics:
        return chemDic({})
    ndics = [nums[ii]*ele for ii, ele in enumerate(dics)]
    return reduce(lambda x, y:x|y, ndics)



def _checkElem(elem: str) -> bool:
    """

    Checks if the string is a atomic element.

    :param elem: param elem: str:
    :param elem: str:
    :param elem: str: 

    """
    if len(elem)>2 or elem[0]!=elem[0].upper() or (len(elem)==2 and elem[1]!=elem[1].lower()):
        return False
    if re.search('[0-9]', elem):
        return False
    if elem.upper() in ["A", "D", "G", "E", "J", "L", "M", "Q", "R", "T", "X", "Z"]:
        return False
    return True


def _checkChemDic(dic: dict) -> bool:
    """

    Checks if the dictionary can be representing a chemical formula.

    :param dic: param dic: dict:
    :param dic: dict:
    :param dic: dict: 

    """
    generics = ("n", "(n+1)", "n-1", "(n-2)", "n-2", "(n-1)", "2n", "4n")
    keys = list(dic.keys())
    vals = list(dic.values())
    checkKeys = False not in [_checkElem(elem) for elem in keys]
    checkVals = False not in [type(ele)==int or ele in generics for ele in vals] 
    return checkKeys and checkVals


def _compDics(dic1: dict, dic2: dict, nn: int=13, dicPoints: dict={}, penalty: float=10) -> float:
    """

    Compares two dictionaries representing chemical formulas.

    :param dic1: dict: 
    :param dic2: dict: 
    :param nn: int:  (Default value = 13)
    :param dicPoints: dict:  (Default value = {})
    :param penalty: float:  (Default value = 10)

    """
    vals1 = 0.#sum([ele if type(ele)==int else _replacenn(ele, nn=nn) for ele in list(dic1.values())]) 
    for elem in dic1:
        poin = 1.
        if elem in dicPoints:
            poin = dicPoints[elem]
        vals1 += poin*dic1[elem]
    #vals2 = sum([ele if type(ele)==int else _replacenn(ele, nn=nn) for ele in list(dic2.values())]) 
    vals2 = 0.
    for elem in dic2:
        poin = 1.
        if elem in dicPoints:
            poin = dicPoints[elem]
        vals2 += poin*dic2[elem]
    tot = vals1+vals2
    l1 = list(dic1)
    l2 = list(dic2)
    pairs = list(set(l1)&set(l2))
    n1 = list(set(l1)-set(l2))
    n2 = list(set(l2)-set(l1))
    points = 0
    for k1 in pairs:
        if k1 in dicPoints:
            poin = dicPoints[k1]
        else:
            poin = 1.
        diff = abs(dic1[k1]-dic2[k1])
        points += abs(diff*poin/tot)
    while n1:
        ele1 = n1.pop()
        iselem = _checkElem(ele1)
        there = False
        if ele1 in dicPoints:
            poin = dicPoints[ele1]
            there = True
        else:
            poin = 1.
        if n2:
            if not iselem: # is not a chemical element
                poins = []
                for ele2 in n2:
                    if not _checkElem(ele2):
                        poins.append( (_compChemNames(ele1, ele2, penalty=penalty)*max(abs(dic1[ele1]-dic2[ele2]), 1.), ele2) )
                if poins:
                    mii = min([ele[0] for ele in poins])
                    im = [ele[0] for ele in poins].index(mii)
                    n2.pop(n2.index(poins[im][1]))
                    points += poin*abs(mii/tot)
                else: # no one to compare
                    if not there:
                        points += dic1[ele1]*abs(penalty/tot)
                    else:
                        points += dic1[ele1]*abs(poin/tot)
            else: # is chemical element
                points += dic1[ele1]*poin/tot
        elif iselem: # no one to compare
            points += dic1[ele1]*poin/tot
        else:
            if not there:
                points += dic1[ele1]*abs(penalty/tot)
            else:
                points += dic1[ele1]*abs(poin/tot)
    for ele2 in n2:
        if _checkElem(ele2):
            if ele2 in dicPoints:
                poin = dicPoints[ele2]
            else:
                poin = 1.
            points += dic2[ele2]*poin/tot
        else:
            if ele2 in dicPoints:
                poin = dicPoints[ele2]
            else:
                poin = penalty
            points += dic2[ele2]*poin/tot            
    return points


def _findBracketPair(stri: str, ini: int=0) -> tuple[int, int]:
    """

    Finds the position of a bracket pair in the string.

    :param stri: str: 
    :param ini: int:  (Default value = 0)

    """
    iret = stri.index("(", ini)
    ii = iret
    jj = stri.index(")", ii+1)
    while "(" in stri[ii+1:jj]:
        ii = stri.index("(", ii+1)
        jj = stri.index(")", jj+1)
    return iret, jj


def _insertElem(elem: str, num: int, dic: dict, ant: str) -> str:
    """

    Inserts an element into a dicionary of chemical formula.

    :param elem: str: 
    :param num: int: 
    :param dic: dict: 
    :param ant: str: 

    """
    if elem=='*':
        there = dic[ant]
        nnum = f'{num}n+{there}'
        dic[ant] = nnum
    elif num:
        if num.isdigit(): # inserting element elem
            dic[elem] = int(num)
        else:
            dic[elem] = num
        ant = elem
    else:
        dic[elem] = 1 # inserting element elem
    return elem


def _readComp(stri: str) -> dict:
    """

    Transforms a string with a compound formula into a dictionary.

    :param stri: str: 

    """
    stri = stri.replace(" ", "")
    elem = stri[0]
    num = ""
    dic = chemDic({})
    ant = ''
    if elem=="(":
        jj, ii = _findBracketPair(stri)
        elem = stri[1:ii]
    else:
        ii = 0
    N = len(stri)
    while ii<N-1:
        ii += 1
        let = stri[ii]
        if let.isdigit():
            if elem!="R":
                num += let
            else:
                elem += let
        elif let=="-" and num:
            num += let
        elif let==let.lower() and let.isalnum():
            if len(elem)==1:
                elem += let
            else:
                num += let
        elif let=="(":
            ant = _insertElem(elem, num, dic, ant)
            kk, mm = _findBracketPair(stri, ini=ii)
            elem = stri[kk+1:mm]
            ii = mm
            num = ""
        elif let==".":
            ant = _insertElem(elem, num, dic, ant)
            if stri[ii+1]=="(":
                kk, mm = _findBracketPair(stri, ini=ii)
                elem = stri[kk+1:mm]
                ii = mm
                num = ""
            else:
                if stri[ii+1].isdigit():
                    num = str(int(stri[ii+1]))
                    step = 1
                else:
                    num = "1"
                    step = 0
                if "." in stri[ii+1:]:
                    jj = stri.index(".", ii+1)
                    elem = stri[ii+1+step:jj]
                    ii = jj-1
                else:
                    elem = stri[ii+1+step:]
                    ii = N
        elif num:
            ant = _insertElem(elem, num, dic, ant)
            elem = let
            num = ""
        else:
            ant = _insertElem(elem, num, dic, ant)
            elem = let
            num = ""
    ant = _insertElem(elem, num, dic, ant)
    return dic





def _readChem(stri: str, magicn: int=13) -> dict:
    """

    Transforms a string with a chemical formula into a dictionary.

    :param stri: str: 
    :param magicn: int:  (Default value = 13)

    """
    stri = stri.replace(" ", "")
    elem = stri[0]
    num = ""
    dic = chemDic({})
    if elem=="(":
        jj, ii = _findBracketPair(stri)
        elem = stri[1:ii]
    else:
        ii = 0
    N = len(stri)
    while ii<N-1:
        ii += 1
        let = stri[ii]
        if let.isdigit():
            if elem!="R":
                num += let
            else:
                elem += let
        elif let=="-" and num:
            num += let
        elif let==let.lower() and let.isalnum():
            if len(elem)==1 and not num:
                elem += let
            else:
                num += let
        elif let=="(":
            if num:
                if num.isdigit(): # inserting element elem
                    dic[elem] = int(num)
                else:
                    dic[elem] = num
            else:
                dic[elem] = 1 # inserting element elem
            kk, mm = _findBracketPair(stri, ini=ii)
            elem = stri[kk+1:mm]
            ii = mm
            num = ""
        elif let==".":
            if num:
                if num.isdigit():
                    dic[elem] = int(num) # inserting element elem
                else:
                    dic[elem] = num
            else:
                dic[elem] = 1 # inserting element elem
            if stri[ii+1]=="(":
                kk, mm = _findBracketPair(stri, ini=ii)
                elem = stri[kk+1:mm]
                ii = mm
                num = ""
            else:
                if stri[ii+1].isdigit():
                    num = str(int(stri[ii+1]))
                    step = 1
                else:
                    num = "1"
                    step = 0
                if "." in stri[ii+1:]:
                    jj = stri.index(".", ii+1)
                    elem = stri[ii+1+step:jj]
                    ii = jj-1
                else:
                    elem = stri[ii+1+step:]
                    ii = N
        elif num:
            if num.isdigit():
                dic[elem] = int(num) # inserting element elem
            else:
                dic[elem] = num
            elem = let
            num = ""
        else:
            dic[elem] = 1 # inserting element elem
            elem = let
            num = ""
    if num:
        if num.isdigit():
            dic[elem] = int(num) # inserting element elem
        else:
            dic[elem] = num
    else:
        dic[elem] = 1 # inserting element elem
    if len(dic)!=1 or len(list(dic)[0])>2:
        done = False
    else:
        done = True
    while not done:
        done = True
        ndic = chemDic(dic.copy())
        for key in ndic.keys():
            kdic = _readChem(key, magicn=magicn)
            if len(kdic)!=1:
                done = False
                num = dic.pop(key)
                if type(num)==str:
                    num = _replacenn(num, magicn)
                dic = dic | num*kdic
    return dic



def binary(num: int) -> list[int]:
    """
    Transforms an integer into a list of binary digits.

    :type num: int
    
    :returns:
    :rtype: list[int]
    """
    digs = []
    while num:
        dig = num%2
        digs.append(dig)
        num = int(num/2)
    return digs


def getPairsTab(reac: ChemReac, dicPoints: dict|None=None, retPoints: bool=False, turn: bool=True, rec: bool=False, red: str="", penalty: float=10.) -> tuple[list[tuple[int, int]], str] | tuple[list[tuple[int, int]], str, float]:
    """
    Returns links for a chemical reaction, along with a table showing
    the scores of performed comparisons.

    If ``retPoints`` is True, also returns the total reaction score
    (zero indicates a balanced reaction).

    This function may call itself recursively depending on reaction
    structure. Some parameters are intended only for internal recursive use.

    :param reac: The chemical reaction.
    :type reac: ChemReac
    :param dicPoints: Dictionary mapping each chemical element to a score.
                      If None, equal weights are used. Default is None.
    :type dicPoints: dict[str, float] | None
    :param retPoints: If True, returns the total reaction score. Default is False.
    :type retPoints: bool
    :param turn: Internal recursive flag controlling reversed reaction testing. Default is True. 
    :type turn: bool
    :param rec: Internal recursive flag. Default is False.
    :type rec: bool
    :param red: Internal recursive parameter specifying element removal. Default is "".
    :type red: str
    :param penalty: Standard mismatch penalty. Default is 10.0.
    :type penalty: float

    :returns: Tuple containing link pairs, scoring table, and optionally
              total reaction score.
    :rtype: tuple[list[tuple[int, int]], str] |
            tuple[list[tuple[int, int]], str, float]
    """
    if dicPoints is None:
            dicPoints = {}
    nn = reac.nn
    mm = reac.mm
    db = reac.db
    #
    pairs = []
    tocheckS = [(ele, reac.Csubs[ii]) for ii, ele in enumerate(reac.dicsubs)]
    snums = reac.lnums[:]
    tocheckP = [(ele, reac.Cprods[ii]) for ii, ele in enumerate(reac.dicprods)]
    pnums = reac.rnums[:]
    torem = list(set([ele for ele in reac.Csubs if ele in reac.Cprods])) # if the same compound in both sides
    ###########################################
    # if only one substrate and many prods, it is surely beeing broken.
    ###########################################
    if not torem and len(tocheckS)==1 and len(tocheckP)>1: # one return here!!
        snums = [ele if type(ele)==int else _replacenn(ele, nn=nn) for ele in snums]
        pnums = [ele if type(ele)==int else _replacenn(ele, nn=nn) for ele in pnums]
        points = {}
        #
        ids = 0
        dic, nums = tocheckS[ids]
        group = [True for ii in tocheckP]
        dicP = chemDic({})
        n = len(tocheckP)
        striTab = ""
        for ii in range(n):
            pnum = pnums[ii]
            dicP = dicP | pnum*tocheckP[ii][0]
            if not striTab:
                striTab += f" & ${tocheckP[ii][0].toFormula()}$"
            else:
                striTab += f" + ${tocheckP[ii][0].toFormula()}$"
        mii = _compDics(dic, dicP, dicPoints=dicPoints, penalty=penalty)
        points[nums] = mii
        striTab += f" \\\\\n${dic.toFormula()}$ & {mii:6.5f}"
        iis = [ii for ii in range(len(tocheckP))]
        iis.sort(reverse=True)
        for ii in iis:
            nump = reac.Cprods[ii]
            if nums<nump:
                pairs.append( (nums, nump) )
            else:
                pairs.append( (nump, nums) )
        toret = pairs, striTab
        if retPoints:
            toret2 = (*toret, points)
            toret = toret2
        return toret
    ###########################################
    # if 1 to remove appearing three times
    # Xn + Xm = Xn+m
    ###########################################
    if len(torem)==1 and (reac.Csubs.count(torem[0])==2 or reac.Cprods.count(torem[0])==2) and reac.Np==2 and reac.Ns==2:
        i1 = torem[0]
        bla = [ele for ele in reac.Csubs if ele!=i1]
        nothingElse=False
        if bla:
            i2 = bla[0]
            insubs = True
        else:
            bla = [ele for ele in reac.Cprods if ele!=i1]
            if bla:
                i2 = bla[0]
                insubs = False
            else:
                nothingElse = True
        if nothingElse:
            point = 0.0
            pairs = [(i1, i1)]
            points = {i1:point}
            striTab = f" & ${reac.getComp(f'{i1}').getChemDic().toFormula()}$ \\\\\n${reac.getComp(f'{i1}').getChemDic().toFormula()}$ & {point:6.5f} \\\\"
            if retPoints:
                return pairs, striTab, points
            else:
                return pairs, striTab
        if i1<i2:
            pairs = [(i1, i2)]
        else:
            pairs = [(i2, i1)]
        if insubs:
            point = _compDics(reac.dicsubs[reac.Csubs.index(i2)], reac.dicprods[reac.Cprods.index(i1)])
        else:
            point = _compDics(reac.dicsubs[reac.Csubs.index(i1)], reac.dicprods[reac.Cprods.index(i2)])
        points = {i1:point, i2:point}
        striTab = f" & ${reac.getComp(f'{i1}').getChemDic().toFormula()}$ \\\\\n${reac.getComp(f'{i2}').getChemDic().toFormula()}$ & {point:6.5f} \\\\"
        if len(snums)!=2 or len(pnums)!=2:
            print(f"Strange reaction!!!!! {reac.ID} - {reac.eq}")
        if retPoints:
            return pairs, striTab, points
        else:
            return pairs, striTab
    ###########################################
    # if 2 to remove, connect them. 
    ###########################################
    if len(torem)==2 and not reac.issueRep:
        i1 = torem[0]
        i2 = torem[1]
        if i1<i2:
            pairs = [(i1, i2)]
        else:
            pairs = [(i2, i1)]
        point = _compDics(reac.dicsubs[reac.Csubs.index(torem[0])], reac.dicprods[reac.Cprods.index(torem[1])])
        points = {i1:point, i2:point}
        striTab = f" & ${reac.getComp(f'{i1}').getChemDic().toFormula()}$ \\\\\n${reac.getComp(f'{i2}').getChemDic().toFormula()}$ & {point:6.5f} \\\\"
        if len(snums)!=2 or len(pnums)!=2:
            print(f"Strange reaction!!!!! {reac.ID} - {reac.eq}")
        if retPoints:
            return pairs, striTab, points
        else:
            return pairs, striTab
    ###########################################
    # must here first check if removing causes trouble. Check numbers
    ###########################################
    passedHere = False
    if not rec and turn and torem:
        ele = torem.pop() # possible returns here
        if torem:
            print("SHOULD BE EMPTY!!!")
        iS = reac.Csubs.index(ele)
        iP = reac.Cprods.index(ele)
        numS = _replacenn(reac.lnums[iS], nn=nn, mm=mm)
        CS = _replacenn(reac.lCcompl[iS], nn=nn, mm=mm)
        numP = _replacenn(reac.rnums[iP], nn=nn, mm=mm)
        CP = _replacenn(reac.rCcompl[iP], nn=nn, mm=mm)
        #nS = _replacenn(numS, nn=nn)*_replacenn(CS, nn=mm)
        #nP = _replacenn(numP, nn=nn)*_replacenn(CP, nn=mm)
        if (numS==numP and CS>CP) or (numS>numP and CS==CP): #nS>nP:#remove only from products
            nreac = reac.copy()
            if CS!=CP and numS==numP: # remove from complement
                nreac.lCcompl[iS] = int(round((CS-CP)/_replacenn(numS, nn=nn, mm=mm), 1))
                nreac.rCcompl[iP] = 0
            elif CS==CP and numS!=numP: # remove from stechiometry
                nreac.lnums[iS] = int(round((numS-numP)/_replacenn(CS, nn=nn, mm=mm), 1))
                nreac.rnums[iP] = 0
            nreac = ChemReac(nreac.rebuildLine(), reduceCompoundS=ele, db=db)
            pairs, striTab, points =  getPairsTab(nreac, retPoints=True, dicPoints=dicPoints, rec=True, red=ele, penalty=penalty) # recursive call with fixed equation
            for i1, i2 in pairs:
                if i1 in reac.Csubs and i1 not in points and i2 in points:
                    points[i1] = points[i2]
                elif i1 in reac.Csubs and i1 in points:
                    #print("Alert?! Check getPairs when Ns>Np and think about it!")
                    pass
                elif i2 in reac.Csubs and i2 not in points and i1 in points:
                    points[i2] = points[i1]
                elif i2 in reac.Csubs and i2 in points:
                    #print("Alert?! Check getPairs when Ns>Np and think about it!")
                    pass
            #pairs = list(set(pairs))
            if retPoints:
                return pairs, striTab, points
            else:
                return pairs, striTab
        elif (numS==numP and CS<CP) or (numS<numP and CS==CP): # nP>nS:#remove only from substrates
            nreac = reac.copy()
            if CS!=CP and numS==numP: # remove from complement
                nreac.lCcompl[iS] = 0
                nreac.rCcompl[iP] = int(round((CP-CS)/_replacenn(numP, nn=nn, mm=mm), 1))
            elif CS==CP and numS!=numP: # remove from stechiometry
                nreac.lnums[iS] = 0
                nreac.rnums[iP] = int(round((numP-numS)/_replacenn(CP, nn=nn, mm=mm), 1))
            #remove only from substrates
            nreac = ChemReac(nreac.rebuildLine(), reduceCompoundS=ele, db=db)
            pairs, striTab, points =  getPairsTab(nreac, retPoints=True, dicPoints=dicPoints, rec=True, red=ele, penalty=penalty) # recursive call with fixed equation
            for i1, i2 in pairs:
                if i1 in reac.Csubs and i1 not in points and i2 in points:
                    points[i1] = points[i2]
                elif i1 in reac.Csubs and i1 in points:
                    #print("Alert?! Check getPairs when Ns>Np and think about it!")
                    pass
                elif i2 in reac.Csubs and i2 not in points and i1 in points:
                    points[i2] = points[i1]
                elif i2 in reac.Csubs and i2 in points:
                    #print("Alert?! Check getPairs when Ns>Np and think about it!")
                    pass
            if retPoints:
                return pairs, striTab, points
            else:
                return pairs, striTab
        else:
            #remove from both
            nreac = reac.copy()
            nreac.lCcompl[iS] = 0
            nreac.rCcompl[iP] = 0
            nreac.lnums[iS] = 0
            nreac.rnums[iP] = 0
            nreac = ChemReac(nreac.rebuildLine(), reduceCompoundS=ele, db=db)
            pairs, striTab, points =  getPairsTab(nreac, retPoints=True, dicPoints=dicPoints, rec=False, red=ele, penalty=penalty) # recursive call with fixed equation
            #passedHere = True # if True reaction might be reverse without need
            #print("****** Strange alert!!!!")
            passedHere = False
            for i1, i2 in pairs:
                if i1 not in points and i2 in points:
                    points[i1] = points[i2]
                elif i2 not in points and i1 in points:
                    points[i2] = points[i1]
            if retPoints:
                return pairs, striTab, points
            else:
                return pairs, striTab
    else:
        torem = []
    points = {}
    ###########################################
    # here the pairs are searched 
    ###########################################
    snums = [ele if type(ele)==int else _replacenn(ele, nn=nn) for ele in snums]
    pnums = [ele if type(ele)==int else _replacenn(ele, nn=nn) for ele in pnums]
    reverted = False
    # below reaction is reverted if more substrates than products
    if len(tocheckS)>len(tocheckP):
        tocheckS, tocheckP = tocheckP, tocheckS
        snums, pnums = pnums, snums
        reverted = True
    # take out zeroed molecules.
    tocheckS = [ele for ele in tocheckS if sum(ele[0].values())]
    tocheckP = [ele for ele in tocheckP if sum(ele[0].values())]
    carbs = [ele[0]["C"] if "C" in ele[0] else 0 for ele in tocheckS]
    carbsP = [ele[0]["C"] if "C" in ele[0] else 0 for ele in tocheckP]
    iis = sorted([ii for ii in range(len(tocheckS))], key=lambda x:sum(tocheckS[x][0].values()), reverse=True)
    #iis = sorted([ii for ii in range(len(tocheckS))], key=lambda x:carbs[x], reverse=True)
    tocheckS = [tocheckS[ii] for ii in iis] # sorts by number of carbons
    snums = [snums[ii] for ii in iis]
    striTab = ""
    while tocheckS and tocheckP:
        resu = []
        first = True
        for ids, (dic, nums) in enumerate(tocheckS):
            snum = snums[ids]
            if snum!=int(snum):
                dicori = snum*dicori
                snum = 1
            else:
                dicori = chemDic(dic.copy())
            for ssnum in range(1, snum+1):
                dic = ssnum*dicori
                ###########################3
                # constructs all combinations of the products
                ###########################3
                n = len(tocheckP)
                groups = [binary(ii) for ii in range(2**n)]
                groups = [ele + (n-len(ele))*[0] for ele in groups]
                groups = [[bool(ele) for ele in ele2] for ele2 in groups]
                groups = [ele for ele in groups if True in ele] 
                comps = []
                allmols = [] # each element in this list is how many molecules of a product are in the comparison
                for igru, group in enumerate(groups):
                    molsused = [[]]
                    dicPs = [chemDic({})]
                    for ii in range(n):
                        if group[ii]:
                            pnum = pnums[ii]
                            #dicP = colapseDics(dicP, multDic(tocheckP[ii][0], pnum))
                            dicPs = [dichere | xx*tocheckP[ii][0] for xx in range(1, pnum+1) for dichere in dicPs]
                            molsused = [ele + [xx] for xx in range(1, pnum+1) for ele in molsused]
                        else:
                            molsused = [ele + [0] for ele in molsused]
                    comps.extend([ ( _compDics(dic, dicP, dicPoints=dicPoints, penalty=penalty), igru, ipn ) for ipn, dicP in enumerate(dicPs) ])
                    allmols += molsused
                ncomps = [ele[0] for ele in comps]
                mii = min(ncomps)
                if first:
                    firstL = " & "
                if ssnum==1:
                    secondL = f"${dicori.toFormula()}$ & "
                else:
                    secondL = f"{ssnum} ${dicori.toFormula()}$ & "
                for ii, mols in enumerate(allmols):
                    if first:
                        plus = False
                        for jj, mol in enumerate(mols):
                            if mol and plus:
                                firstL += " + "
                            elif mol:
                                plus = True
                            if mol>1:
                                firstL += f"{mol} ${tocheckP[jj][0].toFormula()}$"
                            elif mol:
                                firstL += f"${tocheckP[jj][0].toFormula()}$"
                    secondL += f"{ncomps[ii]:6.5f}"
                    if ii!=len(allmols)-1:
                        if first:
                            firstL += " & "
                        secondL += " & "
                if first:
                    striTab += f"{firstL} \\\\\n{secondL}\\\\\n"
                else:
                    striTab += f"{secondL}\\\\\n"
                first = False
                if ncomps.count(mii)>1:
                    iis = [ii for ii, ele in enumerate(ncomps) if ele==mii]
                    tots = [sum(groups[comps[ii][1]]) for ii in iis] # will select the one with the most products
                    imi = iis[tots.index(max(tots))]
                else:
                    imi = ncomps.index(mii)
                mii, igru, ipn = comps[imi]
                resu.append( (mii, groups[igru][:], ipn, allmols[imi], ids, ssnum) )
                #
        ##########################
        # selects the best reconstruction
        ##########################
        comps = [ele[0] for ele in resu]
        mii = min(comps)
        if comps.count(mii)==1:
            imi = comps.index(mii)
        else:
            ismi = [ii for ii in range(len(comps)) if comps[ii]==mii]
            #vals = [sum(multDic(tocheckS[resu[ii][4]][0], resu[ii][5]).values()) for ii in ismi if 0 not in [sum(tocheckP[jj][1].values()) for jj in range(n) if resu[ii][1][jj]]]
            vals = [sum((resu[ii][5]*tocheckS[resu[ii][4]][0]).values()) for ii in ismi]
            imi = ismi[vals.index(max(vals))]
        ##########################
        # removes from substrates
        ##########################
        mii, group, ipn, numsused, ids, ssnum = resu[imi]
        # must check if only one or more molecules reconstructed
        if ssnum==snums[ids]:
            dic, nums = tocheckS.pop(ids)
            bla = snums.pop(ids)
        else:
            dic, nums = tocheckS[ids]
            snums[ids] -= ssnum
        points[nums] = mii
        ##########################
        # removes from products
        ##########################
        iis = [ii for ii in range(n) if group[ii]]
        iis.sort(reverse=True)
        for ii in iis:
            numused = numsused[ii]
            dic2, nump = tocheckP[ii]
            if numused==pnums[ii]:
                dic2, nump = tocheckP.pop(ii)
                bla = pnums.pop(ii)
            else:
                pnums[ii] -= numused
            if nums<nump:
                pairs.append( (nums, nump) )
            else:
                pairs.append( (nump, nums) )
    if (reac.Ns==reac.Np or passedHere) and turn and not rec:
        nreac = reac.reverse(reduceCompoundS=red)
        pairs2, striTab2, points2 = getPairsTab(nreac, retPoints=True, turn=False, dicPoints=dicPoints, penalty=penalty)
        pairs = list(set(pairs)|set(pairs2))
        if striTab2!=striTab:
            striTab += f"\\\\\n{striTab2}"
        for i1, i2 in pairs:
            if i1 in reac.Csubs and i1 not in points and i2 in points2:
                points[i1] = points2[i2]
            elif i1 in reac.Csubs and i1 in points:
                #print("Alert?! Check getPairs when Ns>Np and think about it!")
                pass
            elif i2 in reac.Csubs and i2 not in points and i1 in points2:
                points[i2] = points2[i1]
            elif i2 in reac.Csubs and i2 in points:
                #print("Alert?! Check getPairs when Ns>Np and think about it!")
                pass
        pairs = list(set(pairs))
    toret = pairs, striTab
    if retPoints:
        toret2 = (*toret, points)
        toret = toret2
    return toret



def printReac(rea: KEGGfile|ChemReac, dicPoints: dict|None=None, scaleg: float=0, reverse: bool=False, cstd: int=145) -> str:
    """
    Returns a chemical reaction in LaTeX format as a scheme environment.

    The output string includes:
    - the reaction as a TikZ picture with the links given by the algorthm,
    - scoring tables used during comparison,
    - and the reaction definition itself.

    :param rea: The reaction, either as a KEGGfile object or a ChemReac object.
    :type rea: KEGGfile | ChemReac
    :param dicPoints: Dictionary mapping chemical elements to score values.
                      If None, equal weights are used. Default is None.
    :type dicPoints: dict[str, float] | None
    :param scaleg: Controls scaling of wide tables. If zero, scaling is
                   determined internally via ``cstd``; otherwise a fixed
                   value between 0 and 1 is used. Default is 0.
    :type scaleg: float
    :param reverse: If True, prints the reaction in reversed form.
                    Default is False.
    :type reverse: bool
    :param cstd: Expected number of characters corresponding to page width.
                 Default is 145.
    :type cstd: int

    :returns: A string containing the reaction in LaTeX format.
    :rtype: str
    """
    if scaleg>1.:
        scaleg = 1.
    if dicPoints is None:
            dicPoints = {}
    if type(rea)==KEGGfile:
        if not reverse:
            reac = ChemReac(rea.EQUATION)
        else:
            reac = ChemReac(rea.EQUATION).reverse()
    elif type(rea)==ChemReac:
        if not reverse:
            reac = rea
        else:
            reac = rea.reverse()
    else:
        print("Enter a reaction!")
        return None
    pairs, tab = getPairsTab(reac, dicPoints=dicPoints)
    stri = f"\\bc \\underline{{{rea.ID}}}: ${reac.latex()}$ \\ec"
    kk = [ele for ele in tab.split('\n') if ele] # determine scale here
    nl = len(kk)
    ns = [ele.count('&') for ele in kk]
    nn = ns[0]
    CinCols = (nn+1)*3
    CinNames = [[(len(ele2.strip()), ele2.count("_")) for ele2 in ele.split("&")] for ele in kk if ele.count('&')==nn]
    CinNames = max([sum([ele[0]-ele[1]*2 for ele in ele2]) for ele2 in CinNames])
    if not scaleg:
        scale = min(cstd/(CinNames+CinCols), 1.)
    else:
        scale = scaleg
    stri2 = f'\\bc\n\\scalebox{{{scale}}}{{\\bt{{|c|{(nn)*"c"}|}}\n\\hline\n'
    first = True
    for ii, line in enumerate(kk):
        line = line.replace('\\\\', '')
        nnn = ns[ii]
        if nnn!=nn:
            scale = max(1., float(scaleg))
            nn = nnn
            CinCols = (nn+1)*3
            CinNames = [[(len(ele2.strip()), ele2.count("_")) for ele2 in ele.split("&")] for ele in kk if ele.count('&')==nn]
            CinNames = max([sum([ele[0]-ele[1]*1.6 for ele in ele2]) for ele2 in CinNames])
            if not scaleg:
                scale = min(cstd/(CinNames+CinCols), 1.)
            else:
                scale = scaleg
            stri2 += f'\\hline\n\\et}}\n\\ec\n\\bc\n\\scalebox{{{scale}}}{{\\bt{{|c||{(nn)*"c"}|}}\n\\hline'
            first = True
        stri2 += f"{line} \\\\\n"
        if first:
            stri2 += "\\hline"
            first = False
    stri2 += f'\\hline\n\\et}}\n\\ec'
    stri3 = graphReac(reac, db=reac.db, dicPoints=dicPoints) # try without db
    scale = max(1., float(scaleg))
    striAlmost = f"\\begin{{scheme}}[h]\n\\begin{{framed}}\n{stri}\n{stri2}\n\\bc\n\\scalebox{{{scale}}}{{{stri3}}}\n\n"
    if reac.balanced:
        striBal = "This reaction is balanced! (zero points)"
    else:
        dicL, dicR = reac.dicsLR[0]-reac.dicsLR[1]
        striBal = f"{{\\bf This reaction is not balanced!}} ({reac.points:7.6f} points)\n\n\\bt{{c|c}}\n Missing in Right-HS & Missing in Left-HS \\\\\n\\hline\n $~{dicL.toFormula()}$ &  $~{dicR.toFormula()}$ \\\\\n\\hline\n\\et"
    striFinal = f"{striAlmost}{striBal}\n\\ec\n\\end{{framed}}\n\\caption{{Caption}}\\label{{sc}}\n\\end{{scheme}}"
    return striFinal



##################################
#
# For reading glycan formulas
#
#################################


_mons = [[ele2.strip() for ele2 in ele.split("\t") if ele.strip()] for ele in pkgutil.get_data(__package__, "data/monosachs.txt").decode('utf-8').splitlines() if ele[0]!="#"] # file modified to eliminate Heptose
_dic_mons = {} # dictionary with monosacharides from KEGG
act = ""
for ele in _mons:
    if len(ele)==1:
        act = ele[0]
    elif len(ele)==3:
        _dic_mons[ele[0]] = (ele[1], ele[2], act)
    if len(ele)>3:
        if ele[0]==ele[1]:
            kkk = ele.pop(1)
            _dic_mons[ele[0]] = (ele[1], ele[2], act)
        else:
            _dic_mons[ele[0]] = (ele[2], ele[3], act)
            _dic_mons[ele[1]] = (ele[2], ele[3], act)
            ele[0] = (ele[0], ele.pop(1))


def _removeHifens(ele: str) -> str:
    """

    :param ele: param ele: str:
    :param ele: str:
    :param ele: str: 

    """
    if ele[0]=="-":
        ele = ele[1:]
    if ele and ele[-1]=="-":
        ele = ele[:-1]
    return ele

_monskeys = list(_dic_mons)


def _readGlyc(stri: str) -> list[str]:
    """
    
    Breaks the glycan parts

    :param stri: param stri: str:
    :param stri: str:
    :param stri: str: 

    """
    #if "D/L" in stri:
    #    struc = ["D/L"]
    stri = stri.replace("D/L", "")
    #else:
    struc = []
    was = False
    if "LGalA" in stri:
        was = True
        stri = stri.replace("LGalA", "GalA")
    nums = [(ii, stri.index(ele)) for ii, ele in enumerate(_monskeys) if ele in stri]
    idxs = [ele[1] for ele in nums]
    idxss = list(set([ele[1] for ele in nums]))
    counts = [idxs.count(idd) for idd in idxss]
    for ii, co in enumerate(counts):
        if co>1:
            idx = idxss[ii]
            bla = [ele for ele in nums if ele[1]==idx]
            lengs = [len(_monskeys[ele[0]]) for ele in bla]
            irem = lengs.index(min(lengs))
            rem = bla.pop(irem)
            bla = nums.pop(nums.index(rem))
    theres = [_monskeys[ele[0]] for ele in nums]
    rems = []
    for ii, ele in enumerate(theres):
        for ele2 in [kkk for kkk in theres if kkk!=ele]:
            if ele in ele2:
                rems.append(ii)
    rems = list(set(rems))
    rems.sort(reverse=True)
    bla = [nums.pop(ii) for ii in rems]
    nums = sorted(nums, key=lambda x:x[1])
    ii = 0
    for tup in nums:
        jj = tup[1]
        sac = _monskeys[tup[0]]
        leng = len(sac)
        if jj>ii:
            ele = stri[ii:jj]
            struc.append(ele)
            ii = jj
        if ii+leng>len(stri):
            ele = stri[ii:len(stri)]
            struc.append(ele)
            ii = len(stri)
            break
        struc.append(sac)
        ii += leng
    if len(stri)>ii:
        struc.append(stri[ii:])
    if was:
        ii = struc.index("GalA")
        nstru = ["" for jj in range(len(struc)+1)]
        for jj in range(len(nstru)):
            if jj<ii:
                nstru.append(struc[jj])
            elif jj==ii:
                nstru.append("L")
            else:
                nstru.append(struc[jj-1])
        struc = nstru
    struc = [_removeHifens(ele) for ele in struc if ele]
    return struc


def _numComp(stri: str) -> tuple[str, str|int]:
    """

    Breaks compond and stoichiometry

    :param stri: str: 

    """
    bla = stri.split()
    if bla[0] not in ("n", "(n+1)", "n-1", "(n-2)", "n-2", "(n-1)", "2n", "4n"):
        num = int(bla[0])
    else:
        num = bla[0].replace("(", "").replace(")", "")
    comp = bla[1][:6]
    sub = comp.strip()
    return sub, num



#################################
#
#    stuff
#
#################################



def graphReac(reac: ChemReac, dicPoints: dict={}, colors: list[str]|None=None, penalty: float=10., db: KEGGDB|None=None) -> str:
    """

    Produces the tikz code for showing the chemical reaction with the links produced.

    :param reac: ChemReac: 
    :param dicPoints: dict:  (Default value = {})
    :param colors: list[str]|None:  (Default value = None)
    :param penalty: float:  (Default value = 10.)
    :param db: KEGGDB|None:  (Default value = None)

    """
    if not colors:
        colors = ["red", "blue", "brown", "olive", "orange", "purple", "teal", "violet", "magenta", "cyan", "lime", "green", "pink", "yellow", "black"]
    if type(db)!=KEGGDB:
        db = KEGGDB()
    if db.openned:
        getComp = lambda x: db.allcomps[db.idcomps.index(x)]
    else:
        getComp = lambda x: db.giveMe(x)
    Ns = reac.Ns
    Np = reac.Np
    if Ns<=Np:
        isubs = [getComp(ele) for ele in reac.Csubs]
        ipros = [getComp(ele) for ele in reac.Cprods]
        nsubs = [ii.NAME.split(";")[0] for ii in isubs]
        nprods = [ii.NAME.split(";")[0] for ii in ipros]
        chemS = [ii.getChemDic() for ii in isubs]
        chemP = [ii.getChemDic() for ii in ipros]
        latexS = [reac.latexer(ii) for ii in range(Ns)]
        latexP = [reac.latexer(ii, P=True) for ii in range(Np)]
        lnums = reac.lnums[:]
        rnums = reac.rnums[:]
        for ii, num in enumerate(reac.lCcompl):
            if type(num)!=str:
                lnums[ii] *= num
            else:
                lnums[ii] *= _replacenn(num, nn=reac.nn)
        for ii, num in enumerate(reac.rCcompl):
            if type(num)!=str:
                rnums[ii] *= num
            else:
                rnums[ii] *= _replacenn(num, nn=reac.nn)
        inverted = False
    elif Ns>Np:
        Ns, Np = Np, Ns
        isubs = [getComp(ele) for ele in reac.Cprods]
        ipros = [getComp(ele) for ele in reac.Csubs]
        nsubs = [ii.NAME.split(";")[0] for ii in isubs]
        nprods = [ii.NAME.split(";")[0] for ii in ipros]
        chemS = [ii.getChemDic() for ii in isubs]
        chemP = [ii.getChemDic() for ii in ipros]
        lnums = reac.rnums[:]
        rnums = reac.lnums[:]
        latexS = [reac.latexer(ii, P=True) for ii in range(Ns)]
        latexP = [reac.latexer(ii) for ii in range(Np)]
        for ii, num in enumerate(reac.rCcompl):
            if type(num)!=str:
                lnums[ii] *= num
            else:
                lnums[ii] *= _replacenn(num, nn=reac.nn)
        for ii, num in enumerate(reac.lCcompl):
            if type(num)!=str:
                rnums[ii] *= num
            else:
                rnums[ii] *= _replacenn(num, nn=reac.nn)
        inverted = True
    pairs, tab, points = getPairsTab(reac, dicPoints=dicPoints, retPoints=True, penalty=penalty)
    # problem here if some i in both subs and pros.
    assocPs = [None for ii in range(Np)]
    isubsID = [ele.ID for ele in isubs]
    iprosID = [ele.ID for ele in ipros]
    for i1, i2 in pairs:
        if i1 in isubsID and i2 in iprosID:
            ii = iprosID.index(i2)
            if assocPs[ii]!=None:
                if type(assocPs[ii])==int:
                    assocPs[ii] = [assocPs[ii], isubsID.index(i1)]
                else:
                    assocPs[ii] += [isubsID.index(i1)]
            else:
                assocPs[ii] = isubsID.index(i1)
        elif i1 in iprosID and i2 in isubsID:
            ii = iprosID.index(i1)
            if assocPs[ii]!=None:
                if type(assocPs[ii])==int:
                    assocPs[ii] = [assocPs[ii], isubsID.index(i2)]
                else:
                    assocPs[ii] += [isubsID.index(i2)]
            else:
                assocPs[ii] = isubsID.index(i2)
    assocPs = [ele if ele!=None else -1 for ele in assocPs]
    stri = f"\\node[draw, {colors[0]}] (S1) at (0, 0) {{  ${latexS[0]}$ }};\n"
    for iiss, sub in enumerate(nsubs[1:]):
        iis = iiss+1
        stri += f"\\node (plus{iis}) [right=.1 of S{iis}] {{+}};\n"
        stri += f"\\node[draw, {colors[iis]}] (S{iis+1}) [right=.1 of plus{iis}] {{  ${latexS[iis]}$ }};\n"
    if not inverted:
        stri += f"\\node (A) [right=0.1 of S{Ns}] {{$\\Rightarrow$}};\n"
    else:
        stri += f"\\node (A) [right=0.1 of S{Ns}] {{$\\Leftarrow$}};\n"
    if type(assocPs[0])==list:
        color = "black"
    else:
        color = colors[assocPs[0]]
    stri += f"\\node[draw, {color}] (P1) [right=0.1 of A] {{  ${latexP[0]}$ }};\n"
    for iipp, sub in enumerate(nprods[1:]):
        iip = iipp+1
        if type(assocPs[iip])==list:
            color = "black"
        else:
            color = colors[assocPs[iip]]
        stri += f"\\node (plus{Ns+iip}) [right=.1 of P{iip}] {{+}};\n"
        stri += f"\\node[draw, {color}] (P{iip+1}) [right=.1 of plus{Ns+iip}] {{  ${latexP[iip]}$ }};\n"
    pointsToAna = [ele for ele in points if ele in isubsID]
    flatAssoc = reduce(lambda x,y:x+y, [ele if type(ele)==list else [ele] for ele in assocPs])
    for pointo in pointsToAna:
        ip = isubsID.index(pointo)
        point = points[pointo]
        stri += f"\\node[draw, {colors[ip]}] (V{ip+1}) [below={Ns-ip} of S{ip+1}] {{{point}}};\n"
        stri += f"\\draw[{colors[ip]}] (S{ip+1}) -- (V{ip+1});\n"
        for ii in range(Np):
            if assocPs[ii]==ip or (type(assocPs[ii])==list and ip in assocPs[ii]) or (assocPs[ii]==-1 and ip not in flatAssoc):
                if type(assocPs[ii])==list or assocPs[ii]==-1:
                    color = "black"
                else:
                    color = colors[assocPs[ii]]
                stri += f"\\draw[{color}] (P{ii+1}) -- ($ (P{ii+1}) - (0, {Ns-ip}) $);\n"
                stri += f"\\draw[{color}] ($ (P{ii+1}) - (0, {Ns-ip}) $) -- (V{ip+1});\n"
    return f"\\begin{{tikzpicture}}\n{stri}\\end{{tikzpicture}}"
 





def _guessDic(name: str, magicn: int=13, db: KEGGDB|None=None, forceWeb: bool=False) -> tuple[dict, str]:
    """
    
    Tries to guess what the formula for the compound in str.

    :param name: str: 
    :param magicn: int:  (Default value = 13)
    :param db: KEGGDB|None:  (Default value = None)

    """
    if type(db)!=KEGGDB:
        db = KEGGDB()
    if not forceWeb:
        finder = db.findComp
    else:
        finder = db.findCompWW
    method = ''
    if _checkElem(name):
        dic = _readChem(name)
        solved = True
        method = 'Element'
    elif name=='AH2':
        dic = chemDic({'A': 1, 'H': 2})
        solved = True
        method = 'Known stuff'
    elif name=='R':
        dic = chemDic({'R': 1})
        solved = True
        method = 'Known stuff'
    elif name=='CH3-R':
        dic = chemDic({'C': 1, 'H': 3, 'R': 1})
        solved = True
        method = 'Known stuff'
    elif name=='Sulfur donor':
        dic = chemDic({'S': 1, 'R': 1})
        solved = True
        method = 'Known stuff'
    elif name=='Amino group donor':
        dic = chemDic({'N': 1, 'H':2, 'R': 1})
        solved = True 
        method = 'Known stuff'
    elif name=='ADP-D-ribosyl-acceptor':
        dic = _readChem('(C15H23N5O14P2)n', magicn=magicn) | {'R': 1}
        solved = True
        method = 'Known stuff'
    elif "tRNA" in name:
        if "(" in name:
            kk, mm = _findBracketPair(name)
            brev = name[kk+1:mm]
            if brev in db._abrevs:
                ib = db._abrevs.index(brev)
                dic = _readChem(db._allbrevs[ib][-1])
                dic["tRNA"] = 1
            else:
                dic = chemDic({"tRNA": 1, "R": 1})
        else:
            dic = chemDic({"tRNA": 1, "R": 1})
        if 'thionucleotide' in name:
            dic["S"] = 1
        solved = True
        method = 'tRNA stuff'
    elif "Oxidized" in name or 'Reduced' in name:
        if 'Oxidized' in name:
            term = 'Oxidized'
        else:
            term = 'Reduced'
        key = name.replace(f"{term} ", "")
        comp = finder(key, pprint=False)
        if len(comp)==1:
            comp = db.giveMe(comp[0])
            if "FORMULA" in comp.stuffs:
                dic = _readChem(comp.FORMULA, magicn=magicn)
            else:
                dic = chemDic({key: 1, "e-": 1})
        else:
            dic = chemDic({key: 1, "e-": 1})
        solved = True
        method = 'Redox stuff'
    elif 'electron-transfer' in name.lower():
        nname = name.lower().replace('electron-transfer ', '').strip()
        dic = chemDic({nname: 1})
        method = 'electron-transfer stuff'
        solved = True
    elif 0:#" " not in name:
        comp = db.finder(name, pprint=False, guess=True)
        if comp: comp = db.giveMe(comp[0])
        if comp and "FORMULA" in comp.stuffs:
            dic = _readChem(comp.FORMULA, magicn=magicn)
            solved = True
        else:
            solved = False
    else:
        solved = False
    if not solved: # by guess
        dic = chemDic({name: 1})
        method = 'By name'
    return dic, method







#################################
#
#        Classes
#
#################################





class chemDic(dict):
    """
    Dictionary storing a chemical reaction. Read the manual for details.
    """
    
    
    _generics = ["n", "(n+1)", "n-1", "(n-2)", "n-2", "(n-1)", "2n", "4n"]
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def __rmul__(self, num):
        """Multiplies the values in a dictionary by a number."""
        newdic = chemDic({})
        for key in list(self):
            val = self[key]
            if type(val)!=str:
                nnum = num*self[key]
            else:
                nnum = f"{num}*{val}"
            newdic[key] = nnum
        return newdic
    
    def __or__(self, other):
        """Out of two dictionaries, makes a single one. If the two have a same
        key, the value of the key is summed up."""
        newdic = chemDic({})
        keys = list(self)
        keys2 = list(other)
        for key in keys:
            if key in keys2:
                newdic[key] = self[key] + other[key]
            else:
                newdic[key] = self[key]
        for key in keys2:
            if key not in keys:
                newdic[key] = other[key]
        return newdic
    
    def __sub__(self, other):
        l1 = list(self)
        l2 = list(other)
        pairs = list(set(l1)&set(l2))
        n1 = list(set(l1)-set(l2))
        n2 = list(set(l2)-set(l1))
        dd1 = chemDic({})
        dd2 = chemDic({})
        for ele in n1:
            dd1[ele] = self[ele]
        for ele in n2:
            dd2[ele] = other[ele]
        for ele in pairs:
            v1 = self[ele]
            v2 = other[ele]
            if type(v1)!=str and type(v2)!=str:
                val = v1-v2
                if val>0:
                    dd1[ele] = val
                elif val<0:
                    dd2[ele] = -val
            elif type(v1)!=str:
                val = f"{v2}-{v1}"
                dd2[ele] = val
            elif type(v2)!=str:
                val = f"{v1}-{v2}"
                dd1[ele] = val
            else:
                ii = 109
                while v1==v2:
                    v2 = chr(ii)
                    ii += 1
                val = f"{v1}-{v2}"
                dd1[ele] = val
        return dd1, dd2
    
    def toFormula(self, latex: bool=True) -> str:
        """Transforms dictionary with a chemical formula into a string

        :param latex: Default value = True)
        :param latex: bool:  (Default value = True)

        """
        elems = list(self)
        form = ""
        if "C" in elems:
            num = str(self["C"])
            if num=="1":
                num = ""
            if latex:
                form += "C_{{{:s}}}".format(num)
            else:
                form += "C{:s}".format(num)
            bla = elems.pop(elems.index("C"))
        if "H" in elems:
            num = str(self["H"])
            if num=="1":
                num = ""
            if latex:
                form += "H_{{{:s}}}".format(num)
            else:
                form += "H{:s}".format(num)
            bla = elems.pop(elems.index("H"))
        elems.sort()
        for elem in elems:
            num = str(self[elem])
            if num=="1":
                num = ""
            if latex:
                form += "{:s}_{{{:s}}}".format(elem, num)
            else:
                form += "{:s}{:s}".format(elem, num)
        return form
    
    def check(self) -> bool:
        """ """
        
        keys = list(self.keys())
        vals = list(self.values())
        checkKeys = False not in [_checkElem(elem) for elem in keys]
        checkVals = False not in [type(ele)==int or ele in self._generics for ele in vals] 
        return checkKeys and checkVals




class KEGGDB:
    """
    KEGG Database instance. Please, read the manual for the details.
    """
    
    striKEGG = "http://rest.kegg.jp/%s/%s"
    striKEGGexec = "http://rest.kegg.jp/%s"
    _allbrevs = [("Dol", "Dolichol", "C20H36O(C5H8)n"), ("Cer", "Ceramide", "C19H36NO3R"), ("Ino", "Inositol", "C6H12O6"), ("P", "Phosphoric acid", "H3PO4"), ("S", "Sulfuric acid", "H2SO4"), ("EtN", "Ethanolamine", "C2H7NO"), ("Cho", "Choline", "C5H14NO"), ("Ac", "Acetyl group", "C2H3O"), ("Py", "1-Carboxyethane-1,1-diyl group", "(C3H3O2)2"), ("Me", "Methyl group", "CH3"), ("Gro", "Glycerol", "C3H8O3")] + [("gro", "glycero", "C3H8O3"), ("ery", "erythro", "C4H10O4"), ("thr", "threo", "C4H8O4"), ("rib", "ribo", "C5H10O5"), ("ara", "arabino", "C5H10O5"), ("xyl", "xylo", "C5H10O5"), ("lyx", "lyxo", "C5H10O5"), ("alt", "altro", "C6H12O6"), ("glc", "gluco", "C6H12O6"), ("man", "manno", "C6H12O6"), ("gul", "gulo", "C6H12O6"), ("gal", "galacto", "C6H12O6"), ("tal", "talo", "C6H12O6")] + [("ol", "alditol", "C2H6O2(CH2O)n"), ("onic", "aldonic acid", "C2H4O3(CH2O)n")] + [("methyl", "", "CH3"), ("Fluoro", "Fluoro", "F"), ("Chloro", "Chloro", "Cl"), ("Bromo", "Bromo", "Br"), ("Formyl", "Formyl Group", "CHO"), ("Etn", "Ethanolamine", "C2H7NO"), ("glycer", "Glycerate", "C3H6O4"), ("hydroxy-propanoate", "hydroxy-propanoate", "C3H7O3"), ("Bu", "Butanol", "C4H6O")]
    _dont = [("D", "D-form", ""), ("L", "L-form", ""), ("en", "unsaturated", ""), ("f", "furanose"), ("p", "pyranose"), ("lactone", "lactone"), ("Glc/", "Discarded", "")]
    _removes = [("Anhydro", "", "H2O"), ("d", "deoxy", "O")]
    _aas = [ele.strip().split("\t") for ele in pkgutil.get_data(__package__, "data/AminoAcids.txt").decode('utf-8').splitlines()]
    _aas = [(ele[1], ele[0], ele[4]) for ele in _aas]
    _allbrevs += _aas
    _allbrevs = sorted(_allbrevs, key=lambda x:len(x[0]), reverse=True)
    _abrevs = [ele[0] for ele in _allbrevs]

    def __init__(self, path="", openSTUFF=False, verbose=True):
        try:
            urlopen('https://google.com', timeout=3)
            self.isConnected = True
        except URLError:
            self.isConnected = False
            if verbose: print("Warning: Not connected to the internet. Only local DB available!")
        self.verbose = verbose
        if not path:
            self.pathkegg = str(keggPath)
        else:
            self.pathkegg = path
        self.folders = ["compounds", "glycans", "reactions", "enzymes", 'modules', 'pathways', 'other']
        here = os.listdir(f"{self.pathkegg}")
        for fold in self.folders:
            if fold not in here:
                os.mkdir(f"{self.pathkegg}/{fold}")
        self.comps = [ele for ele in os.listdir(f"{self.pathkegg}/compounds") if re.match(r"C\d{5}", ele)]
        self.glycs = [ele for ele in os.listdir(f"{self.pathkegg}/glycans") if re.match(r"G\d{5}", ele)]
        self.reacs = [ele for ele in os.listdir(f"{self.pathkegg}/reactions") if re.match(r"R\d{5}", ele)]
        self.enzys = [ele for ele in os.listdir(f"{self.pathkegg}/enzymes") if re.match(r"\d\.\d{1,2}\.\d{1,2}\.\d{1,3}", ele)]
        self.modules = [ele for ele in os.listdir(f"{self.pathkegg}/modules") if re.match(r"M\d{5}", ele)]
        self.maps = [ele for ele in os.listdir(f"{self.pathkegg}/pathways") if re.match(r"map\d{5}", ele)]
        #self.comps = os.listdir(f"{self.pathkegg}/compounds")
        #self.glycs = os.listdir(f"{self.pathkegg}/glycans")
        #self.reacs = os.listdir(f"{self.pathkegg}/reactions")
        #self.enzys = os.listdir(f"{self.pathkegg}/enzymes")
        #self.modules = os.listdir(f"{self.pathkegg}/modules")
        #self.maps = os.listdir(f"{self.pathkegg}/pathways")
        self.Ncomps = len(self.comps)
        self.Nglycs = len(self.glycs)
        self.Nreacs = len(self.reacs)
        self.Nenzys = len(self.enzys)
        self.other = os.listdir(f"{self.pathkegg}/other")
        self.openned = False
        if openSTUFF:
            self.redundanceKeys = []
            self.dicRedundances = {}
            self.openSTUFF()
        else:
            self.allreacs = []
            self.idreacs = []
            self.allcomps = []
            self.idcomps = []
            self.dic_comps = {}
            self.compkeys = []
            self.didntmake = []
            self.dicRedundances = {}
            self.redundanceKeys = []
            self.redundanceMain = []
            self.allenzys = []
            self.idenzys = []
            self.orgsInEnz = []
        self.orgsProbReacs = []
        self.orgsProbNets = []
    
    def __repr__(self):
        stri = f"< Comps: {self.Ncomps}, Glycs: {self.Nglycs}, Reacs: {self.Nreacs}, Enzys: {self.Nenzys} >"
        return stri

    def openSTUFF(self):
        """
        Opens lists of compounds, reactions and enzymes as list attributes inside the instance.
        """
        
        verb = self.verbose
        ti = time()
        self.allreacs = [self.giveMe(ele) for ele in self.reacs]
        self.idreacs = [ele.ID for ele in self.allreacs]
        self.allcomps = [self.giveMe(ele) for ele in self.comps+self.glycs]
        self.idcomps = [ele.ID for ele in self.allcomps]
        for ele in self.allreacs:
            ele.web = False
        for ele in self.allcomps:
            ele.web = False
        tf = time()
        if self.verbose: print("Compounds (allcomps, idcomps) and reactions (allreacs, idreacs) read in {:f} s.".format(tf-ti))
        ## enzymes
        ti = tf
        self.allenzys = [self.giveMe(ele) for ele in self.enzys]
        self.idenzys = [ele.ID for ele in self.allenzys]
        self.orgsInEnz = [ele.getOrgs() for ele in self.allenzys]
        if self.orgsInEnz:
            self.orgsInEnzs = list(reduce(lambda x,y:set(x)|set(y), self.orgsInEnz))
        tf = time()
        if self.verbose: print("Enzymes (allenzys, idenzys, orgsInEnz) read in {:f} s.".format(tf-ti))
        ### other stuff
        ti = tf
        self.dic_comps = {} # Dictionary with compound names
        self.compkeys = []
        self.didntmake = []
        for comp in self.allcomps:
            names = [ele.strip() for ele in comp.NAME.lower().split(";")]
            for name in names:
                if name not in self.compkeys:
                    self.dic_comps[name] = comp.ID
                    self.compkeys.append(name)
                else:
                    self.didntmake.append(comp.ID)
        self.openned = True
        tf = time()
        if self.verbose: print("Preparation for findCompsDB in {:f} s.".format(tf-ti))
        ti = tf
        withReps = [ele for ele in self.allcomps if 'REMARK' in ele.stuffs and 'same as:' in ele.REMARK.lower() and not re.search(r'D\d{5}', ele.REMARK)]
        withReps.sort(key=lambda x:x.ID)
        dicSubs = {}
        dones = []
        keys = []
        for ele in withReps:
            lista = re.findall(r'[CG]\d{5}', ele.REMARK)
            lista.append(ele.ID)
            lista.sort()
            comps = [self.giveMe(comp) for comp in lista]
            dics = [comp.getChemDic(nn=13, forceEval=True) for comp in comps]
            ele1 = dics[0]
            nlist = [ele1]
            for ele2 in dics[1:]:
                bla = [ele2!=ele3 for ele3 in nlist]
                if True in bla:
                    nlist.append(ele2)
            if len(nlist)==1 and True not in [ele2 in keys for ele2 in lista]:
                ch = ele.ID
            elif len(dics[0])!=1:
                ch = lista[0]
            else:
                #print(f'{ele} - {nlist}')
                vals = [sum(ele2.values()) for ele2 in dics]
                imax = vals.index(max(vals))
                ch = comps[imax].ID
            appe = False
            for ele2 in lista:
                if ele2 not in keys and ele2!=ch:
                    dicSubs[ele2] = ch
                    keys.append(ele2)
                    appe = True
            if appe and ch not in dones:
                dones.append(ch)
        self.dicRedundances = dicSubs
        self.redundanceKeys = keys
        self.redundanceMain = dones
        tf = time()
        self.verbose = verb
        if self.verbose: print("Searching for redundances took {:f} s.".format(tf-ti))

    def findCompWW(self, stri, onlyequal=True, reteq=False, guess=False, pprint=''):
        """Searches KEGG database in the web for the term.

        :param stri: term searched
        :param onlyequal: If only exact same term should be returned. (Default value = True)
        :param reteq: If should return only the string or the list. (Default value = False)
        :param guess: Should try to guess? (Default value = False)

        """
        if type(pprint)==bool:
            verb = self.verbose
            self.verbose = pprint
        if '/' in stri:
            if self.verbose: print("Warning: Search string should not contain '/'!")
            stri = stri.replace('/', '_')
        stri = stri.lower()
        if ' ' in stri:
            stri.replace(' ', '%20')
        bla = self.find('compound', f'{stri}')
        N = len(bla)
        CS = [ele[0] for ele in bla]
        names = [[ele2.strip().lower() for ele2 in ele[1].split(';')] for ele in bla]
        isthere = [stri in ele for ele in names]
        if self.verbose:
            for ii in range(N):
                print(f"{CS[ii]} -> {'; '.join(names[ii])}")
        toret = CS[:]
        if guess and toret:
            possbs = [self.giveMe(ele) for ele in CS]
            possbs = [ele for ele in possbs if 'generic compound' not in ele.COMMENT.lower() and 'FORMULA' in ele.stuffs]
            forms = [ele.FORMULA for ele in possbs]
            names = [[ele.strip().lower() for ele in ele2.NAME.split(';')] for ele2 in possbs]
            isthere = [stri in ele for ele in names]
            if len(list(set(forms)))==1:
                CS = [possbs[0].ID]
                isthere = [isthere[0]]
            elif names:
                bla = [min([_compChemNames(nstri, nam) for nam in nams]) for nams in names]
                mostSimilar = possbs[bla.index(min(bla))]
                CS = [mostSimilar.ID]
                isthere = [stri in [ele.strip() for ele in mostSimilar.NAME.split(';')]]
            else:
                #print('Could not guess.')
                pass
        if True in isthere and onlyequal:
            toret = [CS[isthere.index(True)]]
        if reteq:
            if True in isthere:
                toret = CS[isthere.index(True)]
            else:
                toret = None
        if type(pprint)==bool:
            self.verbose = verb
        return toret

    def findCompDB(self, stri, ands=True, onlyequal=False, reteq=False, guess=False, pprint=''):
        """Searches the local DB for the term.

        :param stri: term searched
        :param ands: If a list of terms is searched, if all terms should be found. (Default value = True)
        :param onlyequal: If only exact same term should be returned. (Default value = True)
        :param reteq: If should return only the string or the list. (Default value = False)
        :param guess: Should try to guess? (Default value = False)

        """
        if type(pprint)==bool:
            verb = self.verbose
            self.verbose = pprint
        lista = []
        equal = False
        if type(stri)==list:
            nstri = [ele.lower().replace(" ", "") for ele in stri]
            for name in self.dic_comps:
                bla = [ele in name for ele in nstri]
                if ands:
                    checks = False not in bla
                else:
                    checks = True in bla
                if checks:
                    if self.verbose: print("{:s} -> {:s}".format(self.dic_comps[name], name))
                    lista.append(self.dic_comps[name])
        else:
            nstri = stri.lower()
            for name in self.dic_comps:
                if nstri in name:
                    if self.verbose: print("{:s} -> {:s}".format(self.dic_comps[name], name))
                    lista.append(self.dic_comps[name])
                    if nstri==name:
                        equal = True
                        same = self.dic_comps[name]
        if guess and lista and type(stri)!=list:
            possbs = [self.giveMe(ele) for ele in lista]
            possbs = [ele for ele in possbs if 'generic compound' not in ele.COMMENT.lower() and 'FORMULA' in ele.stuffs]
            forms = [ele.FORMULA for ele in possbs]
            #names = [[ele.strip().lower() for ele in ele2.NAME.split(';') if ele.strip().lower()!=nstri] for ele2 in possbs]
            names = [[ele.strip().lower() for ele in ele2.NAME.split(';')] for ele2 in possbs]
            names = [[ele for ele in ele2 if ele!=nstri] for ele2 in names]
            isthere = [nstri in ele for ele in names]
            if len(list(set(forms)))==1:
                lista = [possbs[0].ID]
                equal = True in [isthere[0]]
                same = lista[0]
            elif names:
                bla = [min([_compChemNames(nstri, nam) for nam in nams]) for nams in names]
                mostSimilar = possbs[bla.index(min(bla))]
                lista = [mostSimilar.ID]
                same = lista[0]
                equal = True in [stri in [ele.strip() for ele in mostSimilar.NAME.split(';')]]
            else:
                same = None
                lista = None
        if not lista:
            same = None
        if equal and onlyequal:
            lista = [same]
        if reteq:
            if equal:
                lista = same
            else:
                lista = None
        if type(pprint)==bool:
            self.verbose = verb
        return lista

    def findComp(self, stri, onlyequal=True, reteq=False, guess=False, pprint=''):
        """
        Firts searches the local Database and, if not found, KEGG database via the web.
        """
        l1 = self.findCompDB(stri, onlyequal=onlyequal, reteq=reteq, guess=guess, pprint=pprint)
        if not l1 and self.isConnected:
            l1 = self.findCompWW(stri, onlyequal=onlyequal, reteq=reteq, guess=guess, pprint=pprint)
        #if self.isConnected and web:
        #    l2 = self.findCompWW(self, stri, onlyequal=onlyequal, reteq=reteq, guess=guess, pprint=pprint)
        #else:
        #    l2 = []
        #if not reteq:
        #    if l1 and l2:
        #        l1 = list(set(l1+l2))
        #    elif l2:
        #        l1 = l2
        #else:
        #    if l2 and l1 and l2!=l1:
        #        if self.verbose: print(f"Warning: Found from web ({l2}) and local DB ({l1}) are different. Returning both.")
        #        l1 = [l1, l2]
        return l1

    def giveMe(self, stri):
        """Returns the KEGGfile object.

        :param stri: 

        """
        if re.match(r"C\d{5}", stri):
            if stri in self.comps:
                if self.openned:
                    return self.allcomps[self.idcomps.index(stri)]
                else:
                    return KEGGfile(open(f"{self.pathkegg}/compounds/{stri}").readlines(), db=self)
            else:
                data = self.get(stri)
                print("\n".join(data), file=open(f"{self.pathkegg}/compounds/{stri}", "w"))
                self.Ncomps += 1
                self.comps.append(stri)
                comp = KEGGfile(data, db=self)
                self.allcomps.append(comp)
                self.idcomps.append(stri)
                return comp 
        elif re.match(r"G\d{5}", stri):
            if stri in self.glycs:
                if self.openned:
                    return self.allcomps[self.idcomps.index(stri)]
                else:
                    return KEGGfile(open(f"{self.pathkegg}/glycans/{stri}").readlines(), db=self)
            else:
                data = self.get(stri)
                print("\n".join(data), file=open(f"{self.pathkegg}/glycans/{stri}", "w"))
                self.Nglycs += 1
                self.glycs.append(stri)
                comp = KEGGfile(data, db=self)
                self.allcomps.append(comp)
                self.idcomps.append(stri)
                return comp 
        elif re.match(r"R\d{5}", stri):
            if stri in self.reacs:
                if self.openned:
                    return self.allreacs[self.idreacs.index(stri)]
                else:
                    return KEGGfile(open(f"{self.pathkegg}/reactions/{stri}").readlines(), db=self)
            else:
                data = self.get(stri)
                print("\n".join(data), file=open(f"{self.pathkegg}/reactions/{stri}", "w"))
                self.Nreacs += 1
                self.reacs.append(stri)
                reac = KEGGfile(data, db=self)
                self.allreacs.append(reac)
                self.idreacs.append(stri)
                return KEGGfile(data, db=self)
        elif re.match(r"\d\.\d{1,2}\.\d{1,2}\.\d{1,3}", stri):
            if stri in self.enzys:
                return KEGGfile(open(f"{self.pathkegg}/enzymes/{stri}").readlines(), db=self)
            else:
                data = self.get(stri)
                print("\n".join(data), file=open(f"{self.pathkegg}/enzymes/{stri}", "w"))
                self.Nenzys += 1
                self.enzys.append(stri)
                return KEGGfile(data, db=self)
        elif re.match(r"M\d{5}", stri):
            if stri in self.modules:
                return KEGGfile(open(f"{self.pathkegg}/modules/{stri}").readlines(), db=self)
            else:
                data = self.get(stri)
                print("\n".join(data), file=open(f"{self.pathkegg}/modules/{stri}", "w"))
                return KEGGfile(data, db=self)
        elif re.match(r"map\d{5}", stri):
            if stri in self.maps:
                return KEGGfile(open(f"{self.pathkegg}/pathways/{stri}").readlines(), db=self)
            else:
                data = self.get(stri)
                print("\n".join(data), file=open(f"{self.pathkegg}/pathways/{stri}", "w"))
                return KEGGfile(data, db=self)
        else:
            try:
                data = self.get(stri)
                return KEGGfile(data, db=self)
            except:
                print("Warning: Could not find!!!")
                return None

    def APIexec(self, stri, timeout=30):
        """Executes command directly over the API.

        :param stri: API command to be run. 
        :param timeout: (Default value = 30)

        """
        www = urlopen(self.striKEGGexec%(quote(stri)), timeout=timeout)
        data = www.read().decode()
        return data

    def get(self, stri, timeout=30):
        """Gets the term from KEGG in the web.

        :param stri: term to be gotten.  
        :param timeout: Default value = 30)

        """
        www = urlopen(self.striKEGG%('get', quote(stri)), timeout=timeout)
        data = www.read().decode()
        return data.split('\n')

    def find(self, what, stri):
        """Finds the an object in a set of objects in KEGG database via the web

        :param what: set to the looked 
        :param stri: term to be found there

        """
        www = urlopen(self.striKEGG%('find', quote(f'{what}/{stri}')))
        data = www.read().decode()
        if data.count(':')==data.count('\n'):
            toret = [ele[ele.index(':')+1:].split("\t") for ele in data.split('\n') if ele]
        else:
            toret = [ele.split("\t") for ele in data.split('\n') if ele]
        return toret

    def list(self, what):
        """Lists from the the KEGG database via the web.

        :param what: 

        """
        h = urlopen(self.striKEGG%("list", quote(what)), timeout=30)
        data = h.read()
        data = data.decode("utf-8")
        data = [ele.split("\t") for ele in data.split("\n")]
        data.pop()
        return data
 
    def downloadStuff(self, stuff, forceAll=False, lista=[], logfile=''):
        """
        :param stuff: str: What to download 
        :param forceAll:  bool: donwload new version even if already in database (Default value = False)
        :param lista: list: list of sttring of what to download if none given, downloads it all. (Default value = [])
        :param logfile: str: name (with path) for log file. (Default value = 'logDown.txt')

        """
        if not logfile:
            logfile = 'logDown.txt'
        filLog = open(logfile, 'a')
        fils = os.listdir(self.pathkegg)
        if f"{stuff:s}s" not in fils:
            os.mkdir(f"{self.pathkegg}/{stuff:s}s")
        if not lista:
            ti = time()
            h = urlopen(self.striKEGG%("list", quote(stuff)), timeout=30)
            data = h.read()
            tf = time()
            ttot = tf-ti
            ttot2 = 0
            data = data.decode("utf-8")
            data = [ele.split("\t") for ele in data.split("\n")]
            data.pop()
            bla = [ele[0] for ele in data]
        else:
            bla = lista[:]
            ttot = 0
            ttot2 = 0
        if not forceAll:
            dones = os.listdir(f"{self.pathkegg}/{stuff:s}s")
        else: # enzymes should be downd again for geting new orgs
            dones = []
        todown = [ele for ele in bla if ele not in dones]
        probs = []
        N = len(todown)
        for ii, reac in enumerate(todown):
            ti = time()
            try:
                rea = urlopen(self.striKEGG%("get", "%s"%reac), timeout=30).read()
                print(rea.decode("utf-8"), file=open(f"{self.pathkegg}/{stuff:s}s/{reac:s}", "w"))
                tf = time()
                dt = tf-ti
                ttot += dt
                ttot2 += dt**2
                tbar = ttot/(ii+1)
                std = ttot2 - (ii+1)*tbar**2
                std = (std)**.5/(ii+1)
                left = (N-ii-1)*tbar
                print(f"({ii+1:5d}/{N:5d}) - {reac:s} done ({left/60:5.4f} +- {std/60:5.4f} mins left.)", file=filLog, flush=True)
                if self.verbose: print(f"({ii+1:5d}/{N:5d}) - {reac:s} done in {dt:5.4f} s. ({left/60:5.4f} +- {std/60:5.4f} mins left.)")
            except (HTTPError, URLError) as error:
                tf = time()
                dt = tf-ti
                ttot += dt
                ttot2 += dt**2
                tbar = ttot/(ii+1)
                std = ttot2 - (ii+1)*tbar**2
                std = (std)**.5/(ii+1)
                left = (N-ii-1)*tbar
                probs.append(reac)
                if self.verbose: print(f"** {error:s} ** ({ii+1:5d}/{N:5d}) - {reac:s} done in {dt:5.4fs} s. ({len(probs):4d} probs already) ({left/60:5.4f} +- {std/60:5.4f} mins left.).")
        iprobs = 0
        tttot = ttot
        while probs and iprobs<20:
            iprobs += 1
            todown = probs[:]
            probs = []
            N = len(todown)
            ttot = 0.
            ttot2 = 0.
            for ii, reac in enumerate(todown):
                ti = time()
                try:
                    rea = urlopen(self.striKEGG%("get", "%s"%reac), timeout=30).read()
                    print(rea.decode("utf-8"), file=open(f"{self.pathkegg}/{stuff:s}s/{reac:s}", "w"))
                    tf = time()
                    dt = tf-ti
                    ttot += dt
                    ttot2 += dt**2
                    tbar = ttot/(ii+1)
                    std = ttot2 - (ii+1)*tbar**2
                    std = (std)**.5/(ii+1)
                    left = (N-ii-1)*tbar
                    print(f"({ii+1:5d}/{N:5d}) - {reac:s} done ({left/60:5.4f} mins left.)", file=filLog, flush=True)
                    if self.verbose: print(f"({ii+1:5d}/{N:5d}) - {reac:s} done ({left/60:5.4f} +- {std/60:5.4f} mins left.)")
                except (HTTPError, URLError) as error:
                    tf = time()
                    dt = tf-ti
                    ttot += dt
                    ttot2 += dt**2
                    tbar = ttot/(ii+1)
                    std = ttot2 - (ii+1)*tbar**2
                    std = (std)**.5/(ii+1)
                    left = (N-ii-1)*tbar
                    probs.append(reac)
                    print(f"** {error:s} ** ({ii+1:5d}/{N:5d}) - {reac:s} done in {dt:5.4fs} s. ({len(probs):4d} probs already) ({left/60:5.4f} +- {std/60:5.4f} mins left.).", file=filLog, flush=True)
                    if self.verbose: print(f"** {error:s} ** ({ii+1:5d}/{N:5d}) - {reac:s} done in {dt:5.4fs} s. ({len(probs):4d} probs already) ({left/60:5.4f} +- {std/60:5.4f} mins left.).")
            tttot += ttot
        print(f"Total time donwloading {stuff:s}s: {tttot/60:f} minutes", file=filLog, flush=True)
        filLog.close()
        self.comps = os.listdir(f"{self.pathkegg}/compounds")
        self.Ncomps = len(self.comps)
        self.glycs = os.listdir(f"{self.pathkegg}/glycans")
        self.Nglycs = len(self.glycs)
        self.reacs = os.listdir(f"{self.pathkegg}/reactions")
        self.Nreacs = len(self.reacs)
        self.enzys = os.listdir(f"{self.pathkegg}/enzymes")
        self.Nenzys = len(self.enzys)
        self.modules = os.listdir(f"{self.pathkegg}/modules")
        self.maps = os.listdir(f"{self.pathkegg}/pathways")
        self.other = os.listdir(f"{self.pathkegg}/other")
        if self.openned:
            self.redundanceKeys = []
            self.dicRedundances = {}
            self.openSTUFF()

    def showPair(self, pair, dicPoints={}):
        """
        Prints the pair identification on screen.

        :param pair: 
        :param dicPoints:  (Default value = {})

        """
        comp1 = self.giveMe(pair[0])
        comp2 = self.giveMe(pair[1])
        name1 = comp1.name()
        name2 = comp2.name()
        dic1 = comp1.getChemDic()
        dic2 = comp2.getChemDic()
        points = _compDics(dic1, dic2, dicPoints=dicPoints)
        print(f"{name1:30s} & {name2:30s}\n{dic1} - {dic2} : {points}")

    def showPairs(self, pairs, dicPoints={}):
        """
        Prints a list of pair identifications on screen.

        :param pairs: 
        :param dicPoints:  (Default value = {})

        """
        if pairs:
            for ii, pair in enumerate(pairs):
                comp1 = self.giveMe(pair[0])
                comp2 = self.giveMe(pair[1])
                name1 = comp1.name()
                name2 = comp2.name()
                dic1 = comp1.getChemDic()
                dic2 = comp2.getChemDic()
                points = _compDics(dic1, dic2, dicPoints=dicPoints)
                if pair:
                    print(f"link {ii+1} & {name1} ({pair[0]}) & {name2} ({pair[1]}) & {points}")

    def getReacsForOrg(self, org, forceDown=False):
        """
        Gets list of enzymes and reactions for a given organism.

        :param org: three or four letter code for the organism. 

        """
        if self.isConnected:
            ecs = self.APIexec(f'link/enzyme/{org}').split('\n')
            bla = ecs.pop()
            if bla and self.verbose: print(f"Last new line wasn't empty. Had: < {bla} > !!!")
            ecs = [ele.split('\t')[1].replace('ec:', '') for ele in ecs]
            if forceDown:
                notThere = ecs[:]
            else:
                notThere = [ele for ele in ecs if ele not in self.enzys]
            if notThere: self.downloadStuff('enzyme', lista=notThere)
        else:
            if self.openned:
                ecs = [self.idenzys[ii] for ii in range(self.Nenzys) if org in self.orgsInEnz[ii] ]
            else:
                ecs = [ele for ele in self.enzys if org in self.giveMe(ele).getOrgs()]
        enzs = [self.giveMe(ele) for ele in ecs]
        if enzs:
            ecs = [ele.ID[3:] for ele in enzs]
            # get rns
            rns = reduce(lambda x, y:x+y, [ele.getReacs() for ele in enzs])
            rns = list(set(rns))
            rns.sort()
        else:
            if self.verbose: print(f"Problem: no enzymes found for {org} in local DB. Try to connect to the internet or maybe organism has no enzymes associated to it yet.")
            ecs = []
            rns = []
            self.orgsProbReacs.append(org)
        return ecs, rns
    
    def makeNetForOrg(self, org, remove=[], dicPoints={}, forceDown=False):
        """
        Returns a dictionary with the information and M-matrix concerning the metabolic network for the organism.

        :param org: 
        :param remove:  (Default value = [])
        :param dicPoints: Default value = {})

        """
        pairs = []
        dic_reacPair = {}
        allpoints = {}
        nodes = []
        M = [[]]
        excluded = []
        removeds = []
        redundants = []
        #
        if self.verbose: print('Getting list of enzymes and reactions.')
        ecs, rns = self.getReacsForOrg(org, forceDown=forceDown)
        if self.verbose: print('Getting reactions.')
        if self.isConnected:
            if forceDown:
                notThere = [ele for ele in rns]
            else:
                notThere = [ele for ele in rns if ele not in self.reacs]
            if notThere: self.downloadStuff('reaction', lista=notThere)
        reas = [self.giveMe(ele) for ele in rns]
        reas = [ele for ele in reas if ele]
        dic_data = {}
        if reas:
            reac = None
            chem = None
            if self.verbose: print('Linking compounds via chemical reactions.')
            try:
                #if 1:
                chems = []
                for chem in reas:
                    chems.append( ChemReac(chem.EQUATION, name=chem.ID, db=self, dicPoints=dicPoints) )
                dic_reacPair = {}
                keys = []
                excluded = []
                #allpoints = []
                allpoints = {}
                for reac in chems:
                    pairs, tab, points = getPairsTab(reac, dicPoints=dicPoints, retPoints=True)
                    #allpoints.extend(points.values())
                    for ele in points:
                        if ele in allpoints:
                            allpoints[ele].append(points[ele])
                        else:
                            allpoints[ele] = [points[ele]]
                    nncomps = list(set([ele for pair in pairs for ele in pair]))
                    for pair in pairs:
                        if pair in keys:
                            dic_reacPair[pair].append(reac.ID)
                        else:
                            toexcl = [ele in remove for ele in pair]
                            excl = True in toexcl
                            if not excl:
                                dic_reacPair[pair] = [reac.ID]
                                keys.append(pair)
                            else:
                                excluded.extend([ii for ii in pair if ii in remove])
                #### talvez ainda buscar e excluir: same as, generics, L-, D- stuff
                allcps = reduce(lambda x, y:list(x)+list(y), keys)
                allcps = list(set(allcps))
                allcps.sort()
                redundants = [ele for ele in allcps if ele in self.redundanceKeys]
                pairs = []
                for c1, c2 in keys:
                    change = False
                    c1n = c1
                    c2n = c2
                    if c1 in redundants:
                        change = True
                        c1n = self.dicRedundances[c1]
                    if c2 in redundants:
                        change = True
                        c2n = self.dicRedundances[c2]
                    if change:
                        reacs = dic_reacPair.pop((c1, c2))
                        npair = [c1n, c2n]
                        npair.sort()
                        pair = tuple(npair)
                        dic_reacPair[pair] = reacs
                    else:
                        pair = (c1n, c2n)
                    if pair not in pairs:
                        pairs.append(pair)
                excluded = list(set(excluded))
                allcps = reduce(lambda x, y:list(x)+list(y), pairs)
                allcps = list(set(allcps))
                allcps.sort()
                N = len(allcps)
                dic_nodes = {}
                for ii in range(N):
                    dic_nodes[allcps[ii]] = ii
                #comps = [self.giveMe(comp) for comp in allcps]
                nodes = allcps #[ele.ID for ele in comps]
                K = len(nodes)
                M = [[0 for ii in range(K)] for jj in range(K)]
                nkeys = [list(ele[:]) for ele in pairs]
                for pair in nkeys:
                    ii = nodes.index(pair[0])
                    jj = nodes.index(pair[1])
                    if ii!=jj:
                        M[ii][jj] = 1
                        M[jj][ii] = 1
                # remove zeros
                kis = [sum(ele) for ele in M]
                removeds = []
                torem = list(set([ii for ii in range(K) if kis[ii]==0] + [ii for ii in range(K) if nodes[ii] in excluded]))
                bla = torem.sort(reverse=True)
                for ii in torem:
                    bla = [ele.pop(ii) for ele in M]
                    bla = M.pop(ii)
                    removeds.append(nodes.pop(ii))
                    bla = kis.pop(ii)
                    nodes.pop(ii)
                    K -= 1
                removeds.sort()
            except Exception as error:
                if not reac:
                    reac = chem.ID
                if self.verbose: print(f"Problem! Error ''{error}'' building net for {org} at reaction {reac}!")
                self.orgsProbNets.append( (org, f'Error {error} at reaction {reac}') )
        else:
            self.orgsProbNets.append( (org, 'no reaction') )
            if not pairs: pairs = []
            if not dic_reacPair: dic_reacPair = {}
            if not allpoints: allpoints = {}
            if not nodes: nodes = []
            if not M[0]: M = [[]]
            if not excluded: excluded = []
            if not removeds: removeds = []
            if not redundants: redundants = []
        dic_data = {'ecs': ecs, 'rns':rns, 'pairs': pairs, 'reacByPair': dic_reacPair, 'points': allpoints, 'nodes': nodes, 'Mmatrix': M, 'excluded': excluded, 'removeds': removeds, 'redundants': redundants}
        return dic_data



class KEGGfile:
    """
    Reads a flat-file format of a KEGG object and returns an instance corresponding to it. 
    Read the manual for details.
    """
    
    def __init__(self, data, db=None, forceWeb=False):
        if type(data)==str:
            data = data.split('\n')
        #toprint = False
        struc = "COMMENT"
        other = "No comment"
        stuffs = []
        for line in data:
            if line=="///": # append last data
                if "'" in other and '"' in other:
                    other = other.replace('"', "''")
                if struc not in stuffs:
                    if '"' not in other:
                        exec('self.{:s} = "{:s}"'.format(struc, other))
                    else:
                        exec("self.{:s} = '{:s}'".format(struc, other))
                    stuffs.append(struc)
                elif '"' not in other:
                    if type(eval(f"self.{struc}"))==list:
                        exec(f'self.{struc}.append("{other}")')
                    else:
                        exec(f'self.{struc} = [self.{struc}, "{other}"]')
                else:
                    if type(eval(f"self.{struc}"))==list:
                        exec(f"self.{struc}.append('{other}')")
                    else:
                        exec(f"self.{struc} = [self.{struc}, '{other}']")
                break
            strucn = line[0:12].strip()
            othern = line[12:].strip()
            if strucn and strucn!=struc: # append new data
                if "\\" in other:
                    #print(other)
                    #toprint = True
                    other = other.replace("\\", "\\\\")
                if "'" in other and '"' in other:
                    other = other.replace('"', "''")
                if struc not in stuffs:
                    if '"' not in other:
                        exec(f'self.{struc:s} = "{other:s}"')
                    else:
                        exec(f"self.{struc:s} = '{other:s}'")
                    stuffs.append(struc)
                elif '"' not in other:
                    if type(eval(f"self.{struc}"))==list:
                        exec(f'self.{struc}.append("{other}")')
                    else:
                        exec(f'self.{struc} = [self.{struc}, "{other}"]')
                else:
                    if type(eval(f"self.{struc}"))==list:
                        exec(f"self.{struc}.append('{other}')")
                    else:
                        exec(f"self.{struc} = [self.{struc}, '{other}']")
                struc = strucn
                other = ""
            if other and other[-1]!=";":
                other += "; {:s}".format(othern)
            elif other:
                other += " {:s}".format(othern)
            else:
                other = othern
        ID, TYPE = [ele.strip() for ele in self.ENTRY.split("    ") if ele.strip()]
        #if toprint: print(ID)
        self.ID = ID
        self.TYPE = TYPE
        stuffs.append("ID")
        stuffs.append("TYPE")
        if "NAME" not in stuffs:
            stuffs.append("NAME")
            if "COMPOSITION" in stuffs:
                self.NAME = self.COMPOSITION
            else:
                self.NAME = self.ID
        self.stuffs = stuffs
        if "GENES" in self.stuffs:
            self.orgs = []
        if "ALL_REAC" in self.stuffs:
            self.reacs = []
        if self.TYPE=="Compound" or self.TYPE=="Glycan":
            self.stuffs.append("chemDic")
            self.chemDic = chemDic({})
        if type(db)!=KEGGDB:
            db = KEGGDB()
        if db.openned:
            if not forceWeb:
                self.web = False
            else:
                self.web = True
            self.getComp = lambda x: db.allcomps[db.idcomps.index(x)]
        else:
            self.getComp = lambda x: db.giveMe(x)
            self.web = True
        self.db = db
        if type(self.COMMENT)==list:
            self.COMMENT = self.COMMENT[1]
        self.probs = []
     
    def __repr__(self):
        stri = "< {:s} >".format(self.NAME)
        return stri
    
    def getReacs(self): # for enzymes
        """
        If ALL_REAC is a field of the instance, returns a list of all reactions it contains.
        """
        
        if "ALL_REAC" not in self.stuffs:
            toret = []
        else:
            if self.reacs:
                toret = self.reacs[:]
            else:
                toret = reduce(lambda x, y: x+y, [ele.split() for ele in self.ALL_REAC.split(";")])
                toret = [ele[0:6] for ele in toret if ele[0]=="R"]
        return toret
    
    def getOrgs(self): # for enzymes
        """
        If GENES is a field of the instance, returns all organisms with a gene related to it.
        """
        
        if "GENES" not in self.stuffs:
            toret = []
        else:
            if self.orgs:
                toret = self.orgs[:]
            else:
                toret = [ele[:ele.index(":")].strip().lower() for ele in self.GENES.split(";") if ":" in ele]
        return toret
    
    def name(self):
        """
        Returns a simple string with the name of the object represented by the instance.
        """
        
        self.names = [ele.strip() for ele in self.NAME.split(";")]
        if len(self.names)>1:
            lens = [len(ele.strip()) for ele in self.names]
            return self.names[lens.index(min(lens))]
        else:
            return self.names[0]
     
    def getChemDic(self, nn=13, forceEval=False, viaWeb='bySelf', byRedundance=True):
        """
        If the instance represents a compound of a glycan, tries to reconstruct its chemical formula and returns it in chemDic dictionary form.

        :param nn: Value to substitute generic number n
        :type nn: int (Default value = 13)
        :param forceEval: whether to recompute the dictionary or use the one stored at self.chemDic.
        :type forceEval: bool (Default value = False)
        :param viaWeb: if searches should be done by local database or via the internet.
        :type viaWeb: bool|string (Default value = 'bySelf')
        :param byRedundance: if should look for other compounds/glycans representing the same of the instance.
        :type byRedundance: bool (Default value = True)

        """
        if viaWeb=='bySelf':
            viaWeb=self.web
        if viaWeb:
            findComp = self.db.findCompWW
        else:
            findComp = self.db.findCompDB
        getComp = lambda x: self.db.giveMe(x)
        if self.chemDic and not forceEval:
            return self.chemDic
        self.ByWhat = ''
        if self.TYPE=="Compound":
            if "FORMULA" in self.stuffs:
                dic = _readChem(self.FORMULA, magicn=nn)
                self.hasFormula = True
                self.ByWhat = 'Formula'
            elif self.ID in self.db.redundanceKeys:
                comp = self.db.dicRedundances[self.ID]
                dic = self.db.giveMe(comp).getChemDic(nn=nn, forceEval=forceEval, viaWeb=viaWeb)
                self.ByWhat = f'Redundance - (same as {comp})'
            else:
                dic = chemDic({})
                solved = False
                self.hasFormula = False
                names = self.NAME
                incomments = re.findall(r'C\d{5}', self.COMMENT)
                if incomments:
                    comps = [getComp(ele) for ele in incomments]
                    forms = [_readChem(ele.FORMULA, magicn=nn) for ele in comps if 'FORMULA' in ele.stuffs]
                    if 'dimer' in self.COMMENT and len(forms)==2:
                        dic = _colapseMultipleDics(forms)
                        solved = True
                        self.ByWhat = 'Comments - dimer'
                    elif forms:
                        hasR = [True in [bool(re.search(r'R\d', ele)) or ele=='R' for ele in ele2.keys()] for ele2 in forms]
                        if True in hasR:
                            ii = hasR.index(True)
                            dic = forms[ii]
                            solved = True
                            self.ByWhat = 'Comments - selection (R)'
                        else:
                            nkeys = [len(ele.keys()) for ele in forms]
                            iis = [ii for ii, ele in enumerate(nkeys) if ele==max(nkeys)]
                            if len(iis)==1:
                                dic = forms[iis[0]]
                                solved = True
                                self.ByWhat = 'Comments - selection (most keys)'
                            else:
                                vals = [sum(forms[ii].values()) for ii in iis]
                                ii = vals.index(max(vals))
                                dic = forms[iis[ii]]
                                solved = True
                                self.ByWhat = 'Comments - selection (most atoms)'
                if 'hormone' in names.lower():
                    names = names.split("; ")
                    name = [ele for ele in names if ele.split()[-1].lower()=='hormone']
                    if name:
                        name = name[0]
                        dic[name] = 1
                        solved = True
                        self.ByWhat = 'hormone name'
                    else:
                        names = '; '.join(names)
                if not solved and "Reduced" not in names and 'Oxidized' not in names and not solved:
                    if 'protein' in names:
                        names = names.split('; ')
                        name = [ele for ele in names if 'protein' in ele and ele[len(ele)-7:]=='protein']
                        if name:
                            name = name[0]
                            dic[name] = 1
                            solved = True
                            self.ByWhat = 'protein name'
                        else:
                            names = '; '.join(names)
                knowns = ['tRNA', 'AH2', 'R', 'CH3-R', 'Sulfur donor', 'Amino group donor', 'ADP-D-ribosyl-acceptor']#, 'Electron-transfer']
                there = [ele in self.NAME.split('; ') for ele in knowns]
                if not solved and True in there:
                    names = names.split('; ')
                    ii = there.index(True)
                    dic, method = _guessDic(knowns[ii])
                    self.ByWhat = f'Guess Dic - {method}'
                    solved = True
                if not solved and self.name()=='e-':
                    solved = True
                    dic[self.name()] = 1
                    self.ByWhat = 'Electron'
                if not solved and '-' in self.name() and 'electron-' not in self.name().lower():
                    parts = self.name().split('-')
                    cparts = [self.db.findComp(ele, pprint=False, reteq=True) for ele in parts]
                    nparts = [parts[ii] for ii in range(len(parts)) if not cparts[ii]]
                    cparts = [ele for ele in cparts if ele]
                    for comp in cparts:
                        cp = getComp(comp)
                        dic = dic|cp.getChemDic(nn=nn, forceEval=forceEval, viaWeb=viaWeb)
                    if cparts:
                        for part in nparts:
                            if re.match(r'.* \d{1,2}$', part):
                                newp = re.sub(r'(.*) \d{1,2}$', r'\1', part)
                                comp = self.db.findComp(newp, pprint=False, reteq=True)
                                if comp:
                                    cp = getComp(comp)
                                    dic = dic|cp.getChemDic(nn=nn, forceEval=forceEval, viaWeb=viaWeb)
                            else:
                                dic[part] = 1
                        solved = True
                        self.ByWhat = f'By parts'
                if not solved:
                    dic, method = _guessDic(self.name())
                    self.ByWhat = f'Guess Dic - {method}'
                    solved = True
        elif self.TYPE=='Glycan' and byRedundance and self.ID in self.db.redundanceKeys:
            comp = self.db.dicRedundances[self.ID]
            dic = self.db.giveMe(comp).getChemDic(nn=nn, forceEval=forceEval, viaWeb=viaWeb)
            self.ByWhat = f'Redundance - (same as {comp})'
        elif self.TYPE=="Glycan":
            if "COMPOSITION" in self.stuffs:
                dic = _readComp(self.COMPOSITION)
            else:
                dic = _readComp(self.COMMENT.split(";")[0])
            notFound = []
            found = []
            ndic = chemDic({})
            # first finds the simple composition
            for mon in dic:
                num = dic[mon]
                if type(num)==str:
                    num = _replacenn(num, nn=nn)
                try:
                    comp = self.db.findCompDB(mon, pprint=False, reteq=True)
                except:
                    comp = None
                    self.probs.append(mon)
                if comp:
                    cp = getComp(comp)
                    if 'FORMULA' in cp.stuffs:
                        nndic = num*cp.getChemDic(nn=nn, forceEval=forceEval, viaWeb=viaWeb)  #multDic(getComp(cp).getChemDic(nn=nn), num) # recursive call
                        ndic = ndic | nndic
                        found.append((mon, cp.ID))
            if found:
                for ele, cp in found:
                    num = dic.pop(ele)
            monsThere = []
            monsNotThere = []
            dicThere = {}
            dicNotThere = {}
            for ele in dic: # now will look SNFG
                pmons = _readGlyc(ele)
                num = dic[ele]
                intersec = set(pmons)&set(_monskeys)
                notinter = set(pmons)-set(intersec)
                monsThere.extend(list(intersec))
                monsNotThere.extend([ele for ele in pmons if ele not in intersec])
                for ele2 in intersec:
                    dicThere[ele2] = num
                for ele2 in notinter:
                    dicNotThere[ele2] = num
            for mon in dicThere:
                num = dicThere[mon]
                if type(num)==str:
                    num = _replacenn(num, nn=nn)
                # will try to read chemical formula from each key First from monosacharides
                name = _dic_mons[mon][0]
                try:
                    comp = self.db.findCompDB(mon, pprint=False, reteq=True)
                except:
                    self.probs.append(mon)
                    comp = None
                if comp: cp = getComp(comp)
                if comp and 'FORMULA' in cp.stuffs:
                    nndic = num*cp.getChemDic(nn=nn, forceEval=forceEval, viaWeb=viaWeb) #multDic(getComp(cp).getChemDic(nn=nn), num) # recursive call
                    ndic = ndic | nndic
                    found.append((mon, cp.ID))
                else:
                    listas = []
                    if mon[0]!="-" and mon[-1]!="-":
                        listas = [findComp(stri+name, pprint=False) for stri in ["L-", "D-", ""]]
                    elif mon[0]=="-":
                        listas = [findComp(stri+name[1:], pprint=False) for stri in ["L-", "D-", ""]]
                    elif mon[-1]=="-":
                        listas = [findComp(stri+name[:-1], pprint=False) for stri in ["L-", "D-", ""]]
                    else:
                        listas = [findComp(stri+name[1:-1], pprint=False) for stri in ["L-", "D-", ""]]
                    forms = [[getComp(ele) for ele in lista] for lista in listas]
                    forms = reduce(lambda x, y: x+y, forms)
                    isthere = [name==ele.name() or name in [ele2.strip() for ele2 in ele.NAME.split(';')] or f"D-{name}" in [ele2.strip() for ele2 in ele.NAME.split(';')] or f"L-{name}" in [ele2.strip() for ele2 in ele.NAME.split(';')] for ele in forms]
                    if True in isthere:
                        cp = forms[isthere.index(True)]
                        nndic = num*cp.getChemDic(nn=nn, forceEval=forceEval, viaWeb=viaWeb) #multDic(getComp(cp).getChemDic(nn=nn), num) # recursive call
                        ndic = ndic | nndic
                        found.append((mon, cp.ID))
                    else: # could not find exact term
                        #forms = [ele.FORMULA for ele in forms if 'FORMULA' in ele.stuffs]
                        notFound.append(mon)
                        pass #mores.append(mon)
            iirr = 0
            for mon in dicNotThere: # these were not resolved by SNFG
                num = dicNotThere[mon]
                if type(num)==str:
                    num = _replacenn(num, nn=nn)
                try:
                    pcomp = self.db.findCompDB(mon, pprint=False, reteq=True)
                except:
                    self.probs.append(mon)
                    pcomp = None
                if pcomp:
                    dicb = getComp(pcomp).getChemDic(nn=nn, forceEval=forceEval, viaWeb=viaWeb) #getComp(pcomp).getChemDic() # recursive call
                    nndic = num*dicb
                    ndic = ndic | nndic
                    found.append((mon, pcomp))
                else: # will now try to find chemical formula from abreviations.
                    solved = False
                    monori = mon
                    for abrev in self.db._allbrevs:
                        brev = abrev[0]
                        if brev in mon: # its there!
                            nthere = mon.count(brev)
                            mon = mon.replace(brev, "")
                            chem = abrev[2]
                            dicb = _readChem(chem, magicn=nn)
                            keys = list(dicb)
                            for key in keys:
                                if "," not in key:
                                    dic2 = _readChem(key, magicn=nn)
                                    if len(list(dic2))!=1:
                                        num2 = dicb.pop(key)
                                        if type(num2)!=int:
                                            num2 = _replacenn(num2, nn=nn)
                                        for key2 in list(dic2):
                                            n2 = dic2[key2]
                                            if type(n2)!=int:
                                                n2 = _replacenn(n2, nn=nn)
                                            if key2 in list(dicb):
                                                dicb[key2] += n2*num2
                                            else:
                                                dicb[key2] = n2*num2
                            nndic = num*nthere*dicb
                            ndic = ndic | nndic
                            found.append((monori, brev))
                            solved = True
                    if not solved:
                        monori = mon
                        for abrev in self.db._removes:
                            brev = abrev[0]
                            if brev in mon: # its there!
                                nthere = mon.count(brev)
                                if ',' in mon:
                                    nthere *= mon.count(',')+1
                                mon = mon.replace(brev, "")
                                chem = abrev[2]
                                dicb = _readChem(chem, magicn=nn)
                                keys = list(dicb)
                                for key in keys:
                                    num2 = dicb[key]
                                    if key in list(ndic):
                                        if (ndic[key]-num*nthere*dicb[key])>0:
                                            ndic[key] -= num*nthere*dicb[key]
                                        elif (ndic[key]-num*nthere*dicb[key])==0:
                                            num2 = ndic.pop(key)
                                        else:
                                            pass #probs.append( (dic, key) )
                                found.append((monori, brev))
                                solved = True
                    if not solved and mon:
                        for ele in self.db._dont:
                            mon = mon.replace(ele[0], "")
                        if mon: # some disapear after the replacement here
                            if mon=='*':
                                dic[f"R{iirr}"] = num
                                iirr += 1
                                found.append((mon, "R"))
                            elif mon=="Hep":
                                hepDic = chemDic({"C": 7, "H": 14, "O": 7})
                                ndic = ndic | hepDic
                                found.append((mon, "Heptose"))
                            else:
                                notFound.append(mon)
            self.found = found
            self.notFound = notFound
            self.ByWhat = f'Glycan reconstruction ({len(monsThere)}/{len(monsNotThere)} In/Not - {len(found)} found, {len(notFound)} still unknown).'
            if not ndic:
                dic = {self.NAME: 1}
                self.ByWhat = 'NAME'
            else:
                dic = ndic
        else:
            dic = chemDic({self.NAME: 1})
            self.ByWhat = 'NAME'
        self.chemDic = chemDic(dic.copy())
        return dic
    


class ChemReac:
    """
    Chemical Reaction class. Please read the manual for details.
    """
    
    
    def __init__(self, linha, name="reaction", db=None, nn=13, mm=3, reduceCompoundS='', forceWeb=False, dicPoints={}):
        if type(db)!=KEGGDB:
            db = KEGGDB()
        if db.openned:
            if not forceWeb:
                self.web=False
            else:
                self.web=True
            self.getComp = lambda x: db.allcomps[db.idcomps.index(x)]
        else:
            self.getComp = lambda x: db.giveMe(x)
            self.web = True
        self.db = db
        self.nn = nn
        self.mm = mm
        self.ID = name
        self.eq = linha
        self._reduce = reduceCompoundS
        lhs, rhs = linha.split(" <=> ")
        subs = []
        prods = []
        lnums = []
        rnums = []
        lCcompl = [] # if some compound has a complement beside its name.
        rCcompl = []
        for ele in [ele2.strip() for ele2 in lhs.split(" + ")]:
            if " " not in ele:
                num = 1
                if len(ele)==6:
                    sub = ele
                    compl = 1
                else:
                    sub = ele[0:6]
                    compl = ele[7:-1]
                subs.append(sub)
                lnums.append(num)
                lCcompl.append(compl)
            else:
                sub, num = _numComp(ele)
                if len(sub)==6:
                    compl = 1
                else:
                    compl = sub[7:-1]
                    sub = sub[:6]
                subs.append(sub)
                lnums.append(num)
                lCcompl.append(compl)
        for ele in [ele2.strip() for ele2 in rhs.split(" + ")]:
            if " " not in ele:
                num = 1
                if len(ele)==6:
                    sub = ele
                    compl = 1
                else:
                    sub = ele[0:6]
                    compl = ele[7:-1]
                prods.append(sub)
                rnums.append(num)
                rCcompl.append(compl)
            else:
                sub, num = _numComp(ele)
                if len(sub)==6:
                    compl = 1
                else:
                    compl = sub[7:-1]
                    sub = sub[:6]
                prods.append(sub)
                rnums.append(num)
                rCcompl.append(compl)
        self.Csubs = subs
        self.Cprods = prods
        self.issueRep = False
        # there is the issue of same met twice in subs or prods
        uniqs = list(set(self.Csubs))
        counts = [self.Csubs.count(ele) for ele in uniqs]
        counts = [ii for ii, ele in enumerate(counts) if ele!=1]
        if counts: comp = uniqs[counts[0]]
        must = True
        while must and len(self.Csubs)!=len(uniqs):
            posis = [ii for ii, ele in enumerate(self.Csubs) if ele==comp]
            i1 = posis[0]
            i2 = posis[1]
            cs = lCcompl[i1]
            if cs==lCcompl[i2]:
                self.issueRep = True
                cl = lCcompl.pop(i2)
                self.Csubs.pop(i2)
                nl = lnums.pop(i2)
                ns = lnums[i1]
                if type(nl)!=str and type(ns)!=str:
                    lnums[i1] += nl
                elif type(nl)!=str:
                    lnums[i1] = f'{ns}+{nl}'
                elif type(ns)!=str:
                    lnums[i1] = f'{nl}+{ns}'
                else:
                    ii = 109
                    while ns==nl:
                        nl = chr(ii)
                        ii += 1
                    lnums[i1] = f'{ns}+{nl}'
                if type(cl)!=str and type(cs)!=str:
                    lCcompl[i1] += cl
                elif type(cl)!=str:
                    lCcompl[i1] = f'{cs}+{cl}'
                elif type(cs)!=str:
                    lCcompl[i1] = f'{cl}+{cs}'
                else:
                    ii = 109
                    while cs==cl:
                        cl = chr(ii)
                        ii += 1
                    lCcompl[i1] = f'{cs}+{cl}'
                uniqs = list(set(self.Csubs))
                counts = [self.Csubs.count(ele) for ele in uniqs]
                counts = [ii for ii, ele in enumerate(counts) if ele!=1]
                if counts: comp = uniqs[counts[0]]
            else:
                must = False
        uniqs = list(set(self.Cprods))
        counts = [self.Cprods.count(ele) for ele in uniqs]
        counts = [ii for ii, ele in enumerate(counts) if ele!=1]
        if counts: comp = uniqs[counts[0]]
        must = True
        while must and len(self.Cprods)!=len(uniqs):
            posis = [ii for ii, ele in enumerate(self.Cprods) if ele==comp]
            i1 = posis[0]
            i2 = posis[1]
            cs = rCcompl[i1]
            if cs==rCcompl[i2]:
                self.issueRep = True
                cl = rCcompl.pop(i2)
                self.Cprods.pop(i2)
                nl = rnums.pop(i2)
                ns = rnums[i1]
                if type(nl)!=str and type(ns)!=str:
                    rnums[i1] += nl
                elif type(nl)!=str:
                    rnums[i1] = f'{ns}+{nl}'
                elif type(ns)!=str:
                    rnums[i1] = f'{nl}+{ns}'
                else:
                    ii = 109
                    while ns==nl:
                        nl = chr(ii)
                        ii += 1
                    rnums[i1] = f'{ns}+{nl}'
                if type(cl)!=str and type(cs)!=str:
                    rCcompl[i1] += cl
                elif type(cl)!=str:
                    rCcompl[i1] = f'{cs}+{cl}'
                elif type(cs)!=str:
                    rCcompl[i1] = f'{cl}+{cs}'
                else:
                    ii = 109
                    while cs==cl:
                        cl = chr(ii)
                        ii += 1
                    rCcompl[i1] = f'{cs}+{cl}'
                uniqs = list(set(self.Cprods))
                counts = [self.Cprods.count(ele) for ele in uniqs]
                counts = [ii for ii, ele in enumerate(counts) if ele!=1]
                if counts: comp = uniqs[counts[0]]
            else:
                must = False
        self.isubs = [int(ele[1:]) for ele in self.Csubs] # is it actually used somewhere?
        self.iprods = [int(ele[1:]) for ele in self.Cprods] # is it actually used somewhere?
        self.subsKEGG = [self.getComp(ele) for ele in self.Csubs]
        self.prodsKEGG = [self.getComp(ele) for ele in self.Cprods]
        self.chemS = [ele.getChemDic(viaWeb=self.web) for ele in self.subsKEGG]
        self.chemP = [ele.getChemDic(viaWeb=self.web) for ele in self.prodsKEGG]
        self.lnums = lnums
        self.rnums = rnums
        # below I handle the issue of same met in left and right hand side
        torem = list(set(subs)&set(prods))
        if len(torem)==1:
            cprob = torem[0]
            rem = self.getComp(cprob).getChemDic(viaWeb=self.web)
            isub = subs.index(cprob)
            ipro = prods.index(cprob)
            if lCcompl[isub]==rCcompl[ipro] and ( ('FORMULA' in self.subsKEGG[isub].stuffs and re.search('.*\(.*\)n', self.subsKEGG[isub].FORMULA)) or 'FORMULA' not in self.subsKEGG[isub].stuffs):
                chemS = _colapseMultipleDics([ele.getChemDic(viaWeb=self.web) for jj, ele in enumerate(self.subsKEGG) if jj!=isub], lnums)
                chemP = _colapseMultipleDics([ele.getChemDic(viaWeb=self.web) for jj, ele in enumerate(self.prodsKEGG) if jj!=ipro], rnums)
                totsS = sum([ele for ele in chemS.values()])
                totsP = sum([ele for ele in chemP.values()])
                if len(list(set(chemS)&set(rem)))>len(list(set(chemP)&set(rem))): # in prods should be more
                    lCcompl[isub] = 'n'
                    rCcompl[ipro] = 'n+1'
                elif len(list(set(chemS)&set(rem)))<len(list(set(chemP)&set(rem))): # in subs should be more 
                    lCcompl[isub] = 'n+1'
                    rCcompl[ipro] = 'n'
                elif totsS>totsP:
                    lCcompl[isub] = 'n'
                    rCcompl[ipro] = 'n+1'
                elif totsS<totsP:
                    lCcompl[isub] = 'n+1'
                    rCcompl[ipro] = 'n'
                else: # if arrives here all reactions seem to be ok.
                    pass
        #
        self.lCcompl = lCcompl
        self.rCcompl = rCcompl
        self.allsubs = subs+prods
        self._checks(dicPoints=dicPoints)
        self.Ns = len(self.Csubs)
        self.Np = len(self.Cprods)
    
    def __repr__(self):
        stri = "< {:s} >".format(self.eq)
        return stri
    
    def copy(self, reduceCompoundS=""):
        """
        Returns new instance with the same equation.

        :param reduceCompoundS: Default value = "")

        """
        return ChemReac(self.eq, name=self.ID, nn=self.nn, mm=self.mm, db=self.db, reduceCompoundS=reduceCompoundS)

    def rebuildLine(self):
        """
        Returns string with the reaction stored in the instance.
        """
        
        Ccompl = self.lCcompl
        stri = ""
        for ii in range(self.Ns):
            if ii!=0 and self.lnums[ii] and stri and stri[-1]!=" ":
                stri += " + "
            if self.lnums[ii] and self.lnums[ii]!=1:
                if Ccompl[ii]!=1:
                    stri += f"{self.lnums[ii]} {self.Csubs[ii]}({Ccompl[ii]})"
                else:
                    stri += f"{self.lnums[ii]} {self.Csubs[ii]}"
            elif self.lnums[ii]==1:
                if Ccompl[ii]!=1:
                    stri += f"{self.Csubs[ii]}({Ccompl[ii]})"
                else:
                    stri += f"{self.Csubs[ii]}"
        stri += " <=> "
        Ccompl = self.rCcompl
        for ii in range(self.Np):
            if ii!=0 and self.rnums[ii] and stri[-1:-3]!=" >" and stri[-1]!=" ":
                stri += " + "
            if self.rnums[ii] and self.rnums[ii]!=1:
                if Ccompl[ii]!=1:
                    stri += f"{self.rnums[ii]} {self.Cprods[ii]}({Ccompl[ii]})"
                else:
                    stri += f"{self.rnums[ii]} {self.Cprods[ii]}"
            elif self.rnums[ii]==1:
                if Ccompl[ii]!=1:
                    stri += f"{self.Cprods[ii]}({Ccompl[ii]})"
                else:
                    stri += f"{self.Cprods[ii]}"
        return stri
 
    def reverse(self, reduceCompoundS=""):
        """
        Returns instance with the equation in reverse order.

        :param reduceCompoundS: Default value = "")

        """
        A, B = self.eq.split("<=>")
        return ChemReac(f"{B.strip()} <=> {A.strip()}", name=f"{self.ID} reversed", db=self.db, reduceCompoundS=reduceCompoundS)

    def compNames(self):
        """
        Returns string with the equation where compound codes (CXXXXX) have been substituted by the respective names.
        """
        
        subs = [ele.NAME.split(";")[0] for ele in self.subsKEGG]
        prods = [ele.NAME.split(";")[0] for ele in self.prodsKEGG]
        left = " + ".join(["{:s} {:s}".format(str(self.lnums[ii]), ele) if self.lnums[ii]!=1 else "{:s}".format(ele) for ii, ele in enumerate(subs)])
        right = " + ".join(["{:s} {:s}".format(str(self.rnums[ii]), ele) if self.rnums[ii]!=1 else "{:s}".format(ele) for ii, ele in enumerate(prods)])
        return "{:s} <=> {:s}".format(left, right)

    def latexer(self, ii, P=False):
        """
        Latex form of a compound name.

        :param ii: param P:  (Default value = False)
        :param P: Default value = False)

        """
        if P:
            chems = self.chemP[ii]
            Ccompl = self.rCcompl
            sub = self.prodsKEGG[ii].NAME.split(";")[0]
            numsi = self.rnums
        else:
            chems = self.chemS[ii]
            Ccompl = self.lCcompl
            sub = self.subsKEGG[ii].NAME.split(";")[0]
            numsi = self.lnums
        if numsi[ii]!=1:
            if Ccompl[ii]!=1:
                stri = f"{numsi[ii]} \\underbrace{{(\\textrm{{{sub}$_{{{Ccompl[ii]}}}$}})}}_{{{chems.toFormula()}}}"
            else:
                stri = f"{numsi[ii]} \\underbrace{{\\textrm{{{sub}}}}}_{{{chems.toFormula()}}}"
        else:
            if Ccompl[ii]!=1:
                stri = f"\\underbrace{{\\textrm{{({sub})$_{{{Ccompl[ii]}}}$}}}}_{{{chems.toFormula()}}}"
            else:
                stri = f"\\underbrace{{\\textrm{{{sub}}}}}_{{{chems.toFormula()}}}"
        return stri

    def latex(self):
        """
        Returns a string with equation in latex form.
        """
        
        left = " + ".join([self.latexer(ii) for ii in range(self.Ns)])
        right = " + ".join([self.latexer(ii, P=True) for ii in range(self.Np)])
        return f"{left:s} \\leftrightarrow {right:s}"
    
    def _checks(self, dicPoints={}):
        """

        :param dicPoints: Default value = {})

        """
        nn = self.nn
        mm = self.mm
        if re.findall(r'C\d{5}', self.eq):
            onlyGlyc = False
        else:
            onlyGlyc = True
        subs = [ele.getChemDic(nn=_replacenn(self.lCcompl[jj], nn=nn, mm=mm), forceEval=True, viaWeb=self.web, byRedundance=not onlyGlyc) for jj, ele in enumerate(self.subsKEGG)]
        prods = [ele.getChemDic(nn=_replacenn(self.rCcompl[jj], nn=nn, mm=mm), forceEval=True, viaWeb=self.web, byRedundance=not onlyGlyc) for jj, ele in enumerate(self.prodsKEGG)]
        ########################
        # reduce dictionary of repeated compounds
        ########################
        if self._reduce and self._reduce in self.Csubs and self._reduce in self.Cprods:
            iS = self.Csubs.index(self._reduce)
            iP = self.Cprods.index(self._reduce)
            dicS = subs[iS]
            dicP = prods[iP]
            keys = list(dicS.keys())
            diczero = chemDic(dict(zip(keys, [0 for ii in keys])))
            if sum(dicS.values())>sum(dicP.values()):
                prods[iP] = diczero
                subs[iS] = chemDic(dict(zip(keys, [dicS[key]-dicP[key] for key in keys])))
            else:
                prods[iP] = chemDic(dict(zip(keys, [dicP[key]-dicS[key] for key in keys])))
                subs[iS] = diczero
        self.dicsubs = subs
        self.dicprods = prods
        # put here the complements self.lCcompl, self.rCcompl 
        lnums = [ele if type(ele)!=str else _replacenn(ele, nn=nn) for ele in self.lnums]
        rnums = [ele if type(ele)!=str else _replacenn(ele, nn=nn) for ele in self.rnums]
        sdic = _colapseMultipleDics(subs, lnums)
        pdic = _colapseMultipleDics(prods, rnums)
        self.sdic = sdic
        self.pdic = pdic
        #self.ok = True
        self.points = _compDics(sdic, pdic, dicPoints=dicPoints)
        self.dicsLR = [sdic, pdic]
        if self.points:
            self.balanced = False
        else:
            self.balanced = True


