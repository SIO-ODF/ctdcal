import sys
import db_lib
import db_tables

def main(argv):

    #if the tables are hosed and you need to reset them
    #db_lib.reset_tables()

    #when you start talking to the database
    session = db_lib.start_talk()

    #manual querying of tables
    print(session.query(db_tables.SBE35).all())

    #an example of adding a row to a table one by one or in bulk (from a file, etc)
    #db_lib.setBottleFire(session, 1, 1, 4, 200.5, 100)
    #db_lib.setBottleFireBulk(session, [2,2], [1,1], [1,2], [1020, 900], [20,30])

    #alternate API method of querying whole table
    print(db_lib.getBottleFireAll(session))

    #another row add to table
    #db_lib.setBottleQC(session, 1, 1, 4, 'salt', 3, 'Test comment')
    print(db_lib.getBottleQCAll(session))

    #updating a row in a table
    db_lib.updateBottleQC(session, 1, 1, 4, 'salt', 7)
    print(db_lib.getBottleQCAll(session))

    #when you've finished talking to the database
    db_lib.end_talk(session)

if __name__ == '__main__':
    main(sys.argv[1:])
