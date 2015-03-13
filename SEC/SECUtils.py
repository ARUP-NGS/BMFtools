#!/usr/bin/env python

from lxml import etree
from utilBMF.HTSUtils import ThisIsMadness

"""
Sequencer Error Characterization & Correction (SEC) Utilities.
"""


def ConversionStatsToLaneSuperelement(xmlPath):
    """
    Returns a list of lxml Element objects,
    one for each lane, for all samples.
    Used for error characterization
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


def BuildRunDict(xmlPath):
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
    """
    xmlObj = ConversionStatsToLaneSuperelement(xmlPath)
    newDict = {}
    for lane in xmlObj:
        laneNum = lane.values()[0]
        for tile in lane:
            tileNum = tile.values()[0]
            Raw = tile.getchildren()[0]
            ClusterCount = int(Raw.getchildren()[0].text)
            for entry in range(1, len(Raw)):
                tmpDict = {}
                tmpDict["Yield"] = Raw[i].getchildren()[0].text
                tmpDict["YieldQ30"] = Raw[i].getchildren()[1].text
                tmpDict["QualityScoreSum"] = Raw[i].getchildren()[2].text
                newDict[",".join([laneNum, tileNum, str(entry)])] = tmpDict
    return newDict
