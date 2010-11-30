##
# This script takes an input GFF3 file and produces a useable GBrowse database

import os
import re
import sys
import argparse
import string
import MySQLdb
import subprocess

from shutil import copy
from tempfile import gettempdir
from filelock import FileLock, FileLockException

##
# Regex to find variables in a string so we can replace on them
VARIABLE_RE = re.compile(r'\${([A-Za-z0-9_\.]+)}')

##
# New database stanza wrapper. This stanza goes into the GBrowse master config
gbrowse_db_wrapper = """
[%(SECTION_NAME)s]
description = %(DESCRIPTION)s
path = %(CONF_FILENAME)s"""

class GBrowseConfDuplicateSectionError(Exception):
    """
    A custom exception class that should be raised when attempting to add a duplicate
    section into the GBrowse master conf file
    """
    def __init__(self, value):
        self.error_msg = value
    
    def __str__(self):
        return repr(self.error_msg)        

def buildArgParser():
    """
    Creates an argparser object out of any command line options
    """
    parser = argparse.ArgumentParser(description='Creates a GBrowse database instance ' +
                                                 'from the provided GFF3 file')

    parser.add_argument('-i', '--input_gff3_file', 
                        help='Input GFF3 file that will be loaded into GBrowse')
    parser.add_argument('-s', '--input_sam_file', required=False,
                        help='Optional input SAM file that should be mapped to ' +
                             'the genome loaded into GBrowse')    
    parser.add_argument('-t', '--gbrowse_conf_template', 
                        help='A GBrowse template to generate the data source ' +
                             'configuration from')                        
    parser.add_argument('-o', '--organism', 
                        help='The organism whose data is being loaded into ' +
                             'the GBrowse database')
    parser.add_argument('-g', '--gbrowse_master_conf', default='/opt/gbrowse2/GBrowse.conf',
                        help='Path to GBrowse master configuration file')
    parser.add_argument('-c', '--gbrowse_conf_dir', default='/opt/gbrowse2/',
                        help='Path to GBrowse configuration directory')
    parser.add_argument('-d', '--database',
                        help='Desired name for MySQL database housing annotation data')
    parser.add_argument('-l', '--hostname', help='MySQL server hostname')
    parser.add_argument('-u', '--username', help='Database username')
    parser.add_argument('-p', '--password', help='Database password', default="", required=False)                                  
    args = parser.parse_args()

    return args

def openDatabaseConnection(hostname, username, dbpass):
    """
    Opens a database connection the MySQL server specified using the credentials
    specified.

    Returns a connection to the database in a cursor object
    """
    connection = MySQLdb.connect(host=hostname, user=username, passwd=dbpass)
    return connection

def createAndInitializeDatabase(conn, dbName):
    """
    Utility function that creates a database on the provided host (using he supplied
    username and password) and ensures that permissions are set correctly to allow the 
    gmod and www-data users to work on the database.
    """
    cursor = conn.cursor()
    cursor.execute("CREATE DATABASE %s" % dbName)
    cursor.execute("grant all privileges on %s.* to gmod@localhost" % dbName)
    cursor.execute("grant select on %s.* to 'www-data'@'localhost'" % dbName)
    cursor.close()
 
def createGBrowseDBConfiguration(template, conf_dir, organism, database, conn):                 
    """
    Creates a GBrowse database configuration from the template provided 
    in the input parameters to this script. Uses the organism and database
    parameters to populate some of the fields that must be replaced in the template
    and write a new config file to the GBrowse configuration directory
    """
    rawTemplateStr = open(template).read()
    
    # GBrowse makes use of an inital landmark that the GBrowse page is centered on
    # upon loading an some quick-clickable examples. We want to populate these in our GBrowse
    # configuration by pulling down names from the database
    #
    # TODO: Find a nice way to pull down example landmarks instead of re-using the initial landmark
    initialLandmark = getFeatureByID(conn, database, "1")

    # Replace contents of the template with the parameters passed into this script
    replacedTemplateStr = replaceStr(rawTemplateStr, { "DESCRIPTION": organism, 
                                                       "INITIAL_LANDMARK": initialLandmark[0], 
                                                       "EXAMPLES": initialLandmark[0], 
                                                       "DB_ARGS": "-adaptor DBI::mysql -dsn %s" % database})

    # Now we want to replace the all the necessary variables in the template
    outConf = os.path.join(conf_dir, "%s.conf" % database)
    out = open(outConf, 'w')
    out.write(replacedTemplateStr)
    out.close()    
    
    return outConf
     
def getFeatureByID(conn, db, nameID):
    """
    Pulls down the specified feature from the GBrowse name table

    Returns a string containing the name of the feature
    """
    cursor = conn.cursor()
    cursor.execute("use %s" % db)
    cursor.execute("SELECT name FROM name WHERE id = %s", nameID)
    
    featureName = cursor.fetchone()
    cursor.close()

    return featureName

def replaceStr(str, lookup):
    """Takes a string and replaces all variables in it with values from config"""
    ## Code shamlessly stolen from mmatalka's config.py (vappio library)
    vars = VARIABLE_RE.findall(str)
    if vars:
        for var in vars:
            str = str.replace('${%s}' % var, lookup.get(var))

    return str

def updateGBrowseMasterConfiguration(master_conf, organism, databaseName, dbConf):
    """
    Adds another stanza to the bottom of the GBrowse master configuration
    to ensure that our new database shows up in GBrowse
    
    The stanza uses the following format:
        [<SECTION NAME>]
        description = <description text>
        path = <database configuration filename>
    """
    # We want to check to see if a stanza with the organism name like the one we want to 
    # add to the GBrowse master conf file already exists. If it does we want to raise 
    # an exception of GBrowseConfDuplicateSectionError 
    checkForDuplicateConfSection(master_conf, databaseName)
    
    # Prior to appending to our GBrowse master conf we want to have some insurance
    # in case something screws up so we will copy the conf file to tmp 
    tmpMaster = os.path.join( gettempdir(), os.path.basename(master_conf) )
    copy(master_conf, tmpMaster)
    os.chmod(tmpMaster, 0777) # Need to open up permissions here

    # Implementing a simple form of file locking here to make sure that we 
    # don't run into a race condition while writing to the master conf
    with FileLock(master_conf, timeout=5) as lock:
        master = open(master_conf, 'a')

        # Create new stanza
        master.write( gbrowse_db_wrapper % { 'SECTION_NAME': databaseName,
                                             'DESCRIPTION': organism,
                                             'CONF_FILENAME': os.path.basename(dbConf) } )
        master.close()

def checkForDuplicateConfSection(conf, sectionName):
    """
    Checks a Gbrowse configuration file for the provided section. If this section already
    exists in the configuration files raises a GBrowseConfDuplicateSection exception that 
    is caught in the main function
    """
    confStr = open(conf).read()
        
    if ( confStr.find("[%s]" % sectionName) != -1):
        raise GBrowseConfDuplicateSectionError("%s already exists in the GBrowse master configuration file." % sectionName)           

def cleanUp(conn=None, dbName=None, confFile=None, masterConf=None):
    """
    A cleanup function that deletes the MySQL database or configuration files created
    in the event of a fatal error
    """
    if conn and dbName:
        # Delete the MySQL database if it exists
        cursor = conn.cursor()
        cursor.execute("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = %s", dbName)
        results = cursor.fetchone()
        
        if results:
            cursor.execute("DROP DATABASE %s" % dbName)
            cursor.close()
            conn.close()
    
    if confFile and os.path.isfile(confFile):
        # Delete configuration file if it exists
        os.remove(confFile)
        
    if masterConf:
        tmpMaster = os.path.join(gettempdir(), os.path.basename(masterConf))
        if os.path.isfile(tmpMaster):
            # Restore the original GBrowse master conf file and delete the temp 
            # master conf
            copy(tmpMaster, masterConf)
            os.remove(tmpMaster)

def main(parser):    
    # A valid GFF3 file and its accompanying sequence data (in FASTA) are the 
    # minimum starting point for this script. Additionally a SAM file can be provided
    # that will be mapped against the annotations + sequence loaded.
    dbConn = None
    
    # We want to make sure we have no spaces in our database name
    dbName = parser.database.replace(' ', '_')

    try:
        # Open a connection to the MySQL server for the database operations that will be used
        # throughout this script
        dbConn = openDatabaseConnection(parser.hostname, parser.username, parser.password)

        # Start out by creating the database we will use and setting the proper permissions
        createAndInitializeDatabase(dbConn, dbName)
        
        # Once our DB has been initialized we want to load our GFF3 and FASTA files into the
        # DB
        # TODO: Making a huge assumption here that bioperl and the bp_seqfeature_load script is in 
        # the PATH 
        subprocess.check_call(["bp_seqfeature_load.pl", "-c", "-f", "-d", dbName, 
                                "-u", parser.username,
                                "-p", parser.password,
                                parser.input_gff3_file])

        # Now we need to build our gbrowse config file from the template provided and also append the new
        # entry into the GBrowse master configuration
        confFile = createGBrowseDBConfiguration(parser.gbrowse_conf_template, parser.gbrowse_conf_dir, parser.organism, dbName, dbConn)
        updateGBrowseMasterConfiguration(parser.gbrowse_master_conf, parser.organism, dbName, confFile)
    except MySQLdb.Error, e:
        print "Error initializing database: %s (Error No. %d)" % (e.args[1], e.args[0])
        cleanUp(dbName=dbName, conn=dbConn)
        sys.exit(1)
    except subprocess.CalledProcessError, e:
        print "Error loading GFF3 file into database %s (Error No. %d)" % (parser.database, e.returncode)
        cleanUp(dbName=dbName, conn=dbConn)
        sys.exit(1)
    except IOError, e:
        print "Error writing GBrowse configuration files: %s (Error No. %d)" % (e.args[1], e.args[0])
        cleanUp(dbName=dbName, conn=dbConn, confFile=confFile, masterConf=parser.gbrowse_master_conf)
        sys.exit(1)
    except GBrowseConfDuplicateSectionError, e:
        print "Error adding database to GBrowse master configuration file: %s" % e.error_msg
        cleanUp(dbName=dbName, conn=dbConn, confFile=confFile)
        sys.exit(1)
    except FileLockException:
        print "Error writing to GBrowse master configuration: Could not acquire file lock on file %s" % parser.gbrowse_master_conf
        cleanUp(dbName=parser.database, conn=dbConn, confFile=confFile)
        sys.exit(1)

if __name__ == "__main__":
    main(buildArgParser())

