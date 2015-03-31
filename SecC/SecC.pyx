#!/usr/bin/env python

import logging
from collections import Counter
from operator import attrgetter as oag

import cython
cimport cython
import pysam
from lxml import etree

from utilBMF.HTSUtils import ThisIsMadness, printlog as pl
from MawCluster.PileupUtils import pPileupColumn, PileupReadPair

"""
Sequencer Error Characterization & Correction (SecC) Tools
"""


def ConversionStatsToLaneSuperelement(xmlPath):
    """
    Returns a list of lxml Element objects,
    one for each lane, for all samples.
    To be used for error characterization
    and building an SVM, hopefully.
    """
    return [b for b in etree.parse(xmlPath).getroot()[0].getchildren()
            if b.values()[0] == 'all'][0].getchildren(
                )[0].getchildren()[0]


def GetSequencerInfo(readName):
    """
    Gets read lane, tile, x position, and y position.
    """
    fields = readName.split(":")
    return [int(fields[3]), int(fields[4]), int(fields[5]), int(fields[6])]


def GetRawInfo(xmlObj, runXmlDict=-1, lane=-1, tile=-1):
    """
    Written to take xmlObj as the output of ConversionStatsToLaneElements.
    Returns the Element "Raw" for a given lane and tile.
    """
    if(isinstance(runXmlDict, dict) is False):
        raise ThisIsMadness("a runXmlDict must be provided!")
    if(lane < 0):
        raise ThisIsMadness("Lane must be set to retrieve statistics.")
    if(tile < 0):
        raise ThisIsMadness("Lane must be set to retrieve statistics.")
    laneObj = xmlObj[lane - 1]  # 1 offset for the fact that lanes are 1-based.
    pass


@cython.locals(i=cython.long)
def BuildRunDict(xmlPath, makeGlobal=True):
    """
    Written to take xmlObj as the output of ConversionStatsToLaneElements.
    Returns the a nested dictionary.
    To get a dictionary for a given lane, tile, and read number (NOTE:
    if using a secondary index as an additional read, read 2 is actually
    read 3!), query as follows.
    Yield = newDict[",".join([laneNum, tileNum, str(readNum)])]["Yield"]
    YieldQ30 = newDict[",".join([laneNum, tileNum, str(readNum)])]["YieldQ30"]
    QualityScoreSum = newDict[",".join([laneNum, tileNum,
                                        str(readNum)])]["QualityScoreSum"]
    To access the ClusterCount, use
    ConvXmlDict["laneNum:tileNum:readNum"] with appropriate substitutions
    """
    xmlObj = ConversionStatsToLaneSuperelement(xmlPath)
    if makeGlobal:
        global ConvXmlDict
    ConvXmlDict = {}
    for lane in xmlObj:
        laneNum = lane.values()[0]
        for tile in lane:
            tileNum = tile.values()[0]
            Raw = tile.getchildren()[0]
            ClusterCount = Raw.getchildren()[0].text
            for i in range(1, len(Raw)):
                tmpDict = {}
                tmpDict["Yield"] = Raw[i].getchildren()[0].text
                tmpDict["YieldQ30"] = Raw[i].getchildren()[1].text
                tmpDict["QualityScoreSum"] = Raw[i].getchildren()[2].text
                ConvXmlDict[":".join([laneNum, tileNum, str(i)])] = tmpDict
                ConvXmlDict[":".join(
                    [laneNum, tileNum, str(i), "CC"])] = ClusterCount
    return ConvXmlDict


def GetClusterCount(read, convXml="default"):
    """
    Gets the CC for a read, written for SeqIO.parse
    """
    assert isinstance(convXml, dict)
    rArray = read.name.split(":")
    rNum = read.description.split(" ")[1].split(":")[0]
    return convXml[":".join([rArray[3], rArray[4], rNum, "CC"])]


def GetConvXmlTags(read, convXml="default"):
    """
    Gets the dictionary for those tags for a read, written for SeqIO.parse
    """
    assert isinstance(convXml, dict)
    rArray = read.name.split(":")
    rNum = read.description.split(" ")[1].split(":")[0]
    cx = convXml[":".join([rArray[3], rArray[4], rNum])]
    return {"Y": cx["Yield"], "YQ30": cx["YieldQ30"],
            "QSS": cx["QualityScoreSum"],
            "CC": convXml[":".join([rArray[3], rArray[4], rNum, "CC"])]}