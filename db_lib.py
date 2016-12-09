'''Written in get/set methods for now, should be updated to be more pythonic later.
Need to standardize naming too.

Need uniqueness constraints on rows too

Joseph Gum
December 9, 2016
'''

from sqlalchemy import create_engine
from sqlalchemy.sql import and_
from sqlalchemy.orm import sessionmaker
import db_tables

def start_talk():
    '''Returns a session object that can be used to talk to the db'''
    engine = create_engine('sqlite:///db/odf.db')
    Session = sessionmaker(bind=engine)
    session = Session()
    return session

def end_talk(session):
    session.close()

def reset_tables():
    engine = create_engine('sqlite:///db/odf.db')
    db_tables.Base.metadata.drop_all(engine)
    db_tables.Base.metadata.create_all(engine)
    return True

#Start set/get methods

def setBottleFire(session, in_station, in_cast, in_position, in_pressure, in_time):
    temp = db_tables.BottleFire(station = in_station, cast = in_cast, position = in_position, pressure = in_pressure, time = in_time)
    session.add(temp)
    session.commit()
    return True

def setBottleFireBulk(session, in_station, in_cast, in_position, in_pressure, in_time):
    '''Bulk load instead of add one by one. need to test difference between .add_all() and multiple .adds()
    Assumes inputs other than session are lists to be iterated over separately.
    '''
    temp_list = []
    for a, b, c, d, e in zip(in_station, in_cast, in_position, in_pressure, in_time):
        temp_list.append(db_tables.BottleFire(station = a, cast = b, position = c, pressure = d, time = e))
    for x in temp_list:
        session.add(x)
    session.commit()
    return True

def getBottleFireAll(session):
    return session.query(db_tables.BottleFire).all()

def setSBE35(session, in_station, in_cast, in_position, in_temperature):
    temp = db_tables.SBE35(station = in_station, cast = in_cast, position = in_position, temperature = in_temperature)
    session.add(temp)
    session.commit()
    return True

def setSBE35Bulk(session, in_station, in_cast, in_position, in_temperature):
    '''Bulk load instead of add one by one. need to test difference between .add_all() and multiple .adds()
    Assumes inputs other than session are lists to be iterated over separately.
    '''
    temp_list = []
    for a, b, c, d in zip(in_station, in_cast, in_position, in_temperature):
        temp_list.append(db_tables.SBE35(station = a, cast = b, position = c, temperature = d))
    for x in temp_list:
        session.add(x)
    session.commit()
    return True

def getSBE35All(session):
    return session.query(db_tables.SBE35).all()

def setCoefficientsTemp(session, in_station, in_cast, in_t0, in_t1, in_tp1, in_tp2):
    temp = db_tables.CoefficientsTemp(station = in_station, cast = in_cast, t0 = in_t0, t1 = in_t1, tp1 = in_tp1, tp2 = in_tp2)
    session.add(temp)
    session.commit()
    return True

def setCoefficientsTempBulk(session, in_station, in_cast, in_t0, in_t1, in_tp1, in_tp2):
    '''Bulk load instead of add one by one. need to test difference between .add_all() and multiple .adds()
    Assumes inputs other than session are lists to be iterated over separately.
    '''
    temp_list = []
    for a, b, c, d, e, f in zip(in_station, in_cast, in_t0, in_t1, in_tp1, in_tp2):
        temp_list.append(db_tables.CoefficientsTemp(station = a, cast = b, t0 = c, t1 = d, tp1 = e, tp2 = f))
    for x in temp_list:
        session.add(x)
    session.commit()
    return True

def getCoefficientsTempAll(session):
    return session.query(db_tables.CoefficientsTemp).all()

def setCoefficientsCond(session, in_station, in_cast, in_c0, in_c1, in_c2, in_ct1, in_ct2, in_cp1, in_cp2):
    temp = db_tables.CoefficientsCond(station = in_station, cast = in_cast, c0 = in_c0,
            c1 = in_c1, c2 = in_c2, ct1 = in_ct1, ct2 = in_ct2, cp1 = in_cp1, cp2 = in_cp2)
    session.add(temp)
    session.commit()
    return True

def setCoefficientsCond(session, in_station, in_cast, in_c0, in_c1, in_c2, in_ct1, in_ct2, in_cp1, in_cp2):
    '''Bulk load instead of add one by one. need to test difference between .add_all() and multiple .adds()
    Assumes inputs other than session are lists to be iterated over separately.
    '''
    temp_list = []
    for a, b, c, d, e, f, g, h, i in zip(in_station, in_cast, in_c0, in_c1, in_c2, in_ct1, in_ct2, in_cp1, in_cp2):
        temp_list.append(db_tables.CoefficientsCond(station = a, cast = b, c0 = c, c1 = d, c2 = e, ct1 = f, ct2 = g, cp1 = h, cp2 = i))
    for x in temp_list:
        session.add(x)
    session.commit()
    return True

def getCoefficientsCondAll(session):
    return session.query(db_tables.CoefficientsCond).all()

def setCoefficientsOxygen(session, in_station, in_cast, ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9):
    temp = db_tables.CoefficientsOxygen(station = in_station, cast = in_cast, c1 = ic1, c2 = ic2,
            c3 = ic3, c4 = ic4, c5 = ic5, c6 = ic6, c7 = ic7, c8 = ic8, c9 = ic9)
    session.add(temp)
    session.commit()
    return True

def setCoefficientsOxygen(session, in_station, in_cast, ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9):
    '''Bulk load instead of add one by one. need to test difference between .add_all() and multiple .adds()
    Assumes inputs other than session are lists to be iterated over separately.
    '''
    temp_list = []
    for a, b, c, d, e, f, g, h, i, j, k in zip(in_station, in_cast, ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9):
        temp_list.append(db_tables.CoefficientsOxygen(station = a, cast = b, c1 = c, c2 = d, c3 = e, c4 = f, c5 = g, c6 = h, c7 = i, c8 = j, c9 = k))
    for x in temp_list:
        session.add(x)
    session.commit()
    return True

def getCoefficientsOxygenAll(session):
    return session.query(db_tables.CoefficientsOxygen).all()

def setBottleQC(session, in_station, in_cast, in_position, in_codeproperty, in_qc, in_comment):
    temp = db_tables.BottleQC(station = in_station, cast = in_cast, position = in_position, codeproperty = in_codeproperty, qc = in_qc, comment = in_comment)
    session.add(temp)
    session.commit()
    return True

def setBottleQCBulk(session, in_station, in_cast, in_position, in_codeproperty, in_qc, in_comment):
    '''Bulk load instead of add one by one. need to test difference between .add_all() and multiple .adds()
    Assumes inputs other than session are lists to be iterated over separately.
    '''
    temp_list = []
    for a, b, c, d, e, f in zip(in_station, in_cast, in_position, in_codeproperty, in_qc, in_comment):
        temp_list.append(db_tables.SBE35(station = a, cast = b, position = c, codeproperty = d, qc = e, comment = f))
    for x in temp_list:
        session.add(x)
    session.commit()
    return True

def updateBottleQC(session, in_station, in_cast, in_position, in_codeproperty, in_qc, in_comment=False):

    temp = session.query(db_tables.BottleQC).filter(
        and_(
            db_tables.BottleQC.station == in_station,
            db_tables.BottleQC.cast == in_cast,
            db_tables.BottleQC.position == in_position,
            db_tables.BottleQC.codeproperty == in_codeproperty,
            )
            ).first()
    temp.qc = in_qc
    if in_comment:
        temp.comment = in_comment
    session.commit()
'''    if in_comment == False:
        temp = session.query().\
            filter(db_tables.BottleQC.station == in_station and db_tables.BottleQC.cast == in_cast
                and db_tables.BottleQC.position == in_position and db_tables.BottleQC.codeproperty == in_codeproperty)
        temp.qc = in_qc
    else:
        session.query().\
            filter(db_tables.BottleQC.station == in_station and db_tables.BottleQC.cast == in_cast
                and db_tables.BottleQC.position == in_position and db_tables.BottleQC.codeproperty == in_codeproperty).\
                update({"qc":in_qc, "comment":in_comment})
    session.commit()'''

def getBottleQCAll(session):
    return session.query(db_tables.BottleQC).all()
