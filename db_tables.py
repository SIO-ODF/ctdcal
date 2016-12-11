'''
Needs to be discussed, normalization? better schema design? Tables currently
reflect what is expected in the report post cruise.

Needs to be squared away with bottle_db database design to cutdown on redundant
tables, must be fixed post-cruise
'''

from sqlalchemy import Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

'''class WireDeployment(Base):
    __tablename__ = 'wire_deployment'

    id = Column(Integer, primary_key = True)
    station = Column(Integer)
    cast = Column(Integer)
    deployment_type = Column(String) #not limited, but validated against table of ctd, nettow, etc
'''


class SBE35(Base):
    __tablename__ = 'SBE35'

    id = Column(Integer, primary_key = True)
    station = Column(Integer)
    cast = Column(Integer)
    position = Column(Integer)
    temperature = Column(Float)
    #ForeignKeyConstraint(['station', 'cast', 'position'], ['bottle_fire.station', 'bottle_fire.cast', 'bottle_fire.position'])

    def __repr__(self):
        return "<SBE35(station = '%s', cast = '%s', position = '%s', temperature = '%s')" %(
            self.station, self.cast, self.position, self.temperature)

    def __odfFormat(self):
        return "'%s'/'%s' '%s'" %(self.station, self.cast, self.cast*100 + self.position)


class BottleFire(Base):
    __tablename__ = 'bottle_fire'

    id = Column(Integer, primary_key = True)
    station = Column(Integer)
    cast = Column(Integer)
    position = Column(Integer)
    pressure = Column(Float)
    time = Column(Integer)

    def __repr__(self):
        return "<Bottle(station = '%s', cast = '%s', position = '%s', pressure = '%s', time = '%s')" %(
            self.station, self.cast, self.position, self.pressure, self.time)

class CoefficientsTemp(Base):
    __tablename__ = 'coefficients_temp'

    id = Column(Integer, primary_key = True)
    station = Column(Integer)
    cast = Column(Integer)
    t0 = Column(Float)
    t1 = Column(Float)
    tp1 = Column(Float)
    tp2 = Column(Float)

    def __repr__(self):
        return "<Coefficients Temp(station = '%s', cast = '%s', t0 = '%s', t1 = '%s', tp1 = '%s', tp2 = '%s')" %(
            self.station, self.cast, self.t0, self.t1, self.tp1, self.tp2)

class CoefficientsCond(Base):
    __tablename__ = 'coefficients_cond'

    id = Column(Integer, primary_key = True)
    station = Column(Integer)
    cast = Column(Integer)
    c0 = Column(Float)
    c1 = Column(Float)
    c2 = Column(Float)
    ct1 = Column(Float)
    ct2 = Column(Float)
    cp1 = Column(Float)
    cp2 = Column(Float)

    def __repr__(self):
        return "<Coefficients Cond(station = '%s', cast = '%s', c0 = '%s', c1 = '%s', c2 = '%s', ct1 = '%s', ct2 = '%s', cp1 = '%s', cp2 = '%s')" %(
            self.station, self.cast, self.c0, self.c1, self.c2, self.ct1, self.ct2, self.cp1, self.cp2)

class CoefficientsOxygen(Base):
    __tablename__ = 'coefficients_oxygen'

    id = Column(Integer, primary_key = True)
    station = Column(Integer)
    cast = Column(Integer)
    c1 = Column(Float)
    c2 = Column(Float)
    c3 = Column(Float)
    c4 = Column(Float)
    c5 = Column(Float)
    c6 = Column(Float)
    c7 = Column(Float)
    c8 = Column(Float)
    c9 = Column(Float)

    def __repr__(self):
        return "<Coefficients Cond(station = '%s', cast = '%s', c1 = '%s', c2 = '%s', c3 = '%s', c4 = '%s', c5 = '%s', c6 = '%s', c7 = '%s', c8 = '%s', c9 = '%s')" %(
            self.station, self.cast, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6, self.c7, self.c8, self.c9)

class BottleQC(Base):
    __tablename__ = 'bottle_qc'

    id = Column(Integer, primary_key = True)
    station = Column(Integer)
    cast = Column(Integer)
    position = Column(Integer)
    codeproperty = Column(String)
    qc = Column(Integer)
    comment = Column(String)

    def __repr__(self):
        return "<Bottle QC(station = '%s', cast = '%s', position = '%s', codeproperty = '%s', qc = '%s', comment = '%s')" %(
            self.station, self.cast, self.position, self.codeproperty, self.qc, self.comment)
