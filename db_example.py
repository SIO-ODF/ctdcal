import sys
import db_lib
import db_tables

def main(argv):
    #db_lib.reset_tables()
    session = db_lib.start_talk()
    print(session.query(db_tables.SBE35).all())
    #db_lib.setBottleFire(session, 1, 1, 4, 200.5, 100)
    #db_lib.setBottleFireBulk(session, [2,2], [1,1], [1,2], [1020, 900], [20,30])
    print(db_lib.getBottleFireAll(session))
    #db_lib.setBottleQC(session, 1, 1, 4, 'salt', 3, 'Test comment')
    print(db_lib.getBottleQCAll(session))
    db_lib.updateBottleQC(session, 1, 1, 4, 'salt', 7, 'test comment')
    print(db_lib.getBottleQCAll(session))
    db_lib.end_talk(session)

if __name__ == '__main__':
    main(sys.argv[1:])
