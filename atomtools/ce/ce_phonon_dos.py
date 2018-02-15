import sqlite3 as sq
import numpy as np
import io
from ase.db import connect

class PhononDOS_DB(object):
    def __init__( self, db_name ):
        self.db_name = db_name
        self.check_db()

    def check_db(self):
        """
        Creates the required table
        """
        required_table = "phonon_dos"
        required_fields = ["id","name","atID","omega_e","dos_e"]
        types = {
            "id":"integer",
            "name":"text",
            "atID":"integer",
            "omega_e":"blob",
            "dos_e":"blob"
        }

        conn = sq.connect( self.db_name )
        cur = conn.cursor()

        sql = "create table if not exists %s (id integer)"%(required_table)
        cur.execute(sql)
        conn.commit()

        for field in required_fields:
            try:
                sql = "alter table %s add column %s %s"%(required_table,field,types[field])
                cur.execute(sql)
            except Exception as exc:
                pass
        conn.commit()
        conn.close()

    def get_max_id( self ):
        """
        Returns the current maximum id in the DB
        """
        conn = sq.connect( self.db_name )
        cur = conn.cursor()
        sql = "select id from phonon_dos"
        cur.execute(sql)
        all_entries = cur.fetchall()
        conn.close()
        maxID = -1
        for entry in all_entries:
            if ( entry[0] > maxID ):
                maxID = entry[0]
        return maxID

    def save( self, name=None, atID=None, omega_e=None, dos_e=None ):
        if ( name is None ):
            raise ValueError( "No name specified!" )
        if ( atID is None ):
            raise ValueError( "No atom id specified!" )
        if ( omega_e is None ):
            raise ValueError( "No energies/frequencies specified!" )
        if ( dos_e is None ):
            raise ValueError( "No density of states specified!" )

        omega_blob = self.array_to_blob(omega_e)
        dos_blob = self.array_to_blob(dos_e)
        newID = self.get_max_id()+1

        conn = sq.connect( self.db_name )
        cur = conn.cursor()
        sql = "insert into phonon_dos (id,atID,name,omega_e,dos_e) values (?,?,?,?,?)"
        entries = (newID,atID,name,omega_blob,dos_blob)
        cur.execute(sql,entries)
        conn.commit()
        conn.close()

    def get_all( self ):
        """
        Returns a list of all
        """
        maxID = self.get_max_id()
        allres = []
        for uid in range(maxID):
            try:
                allres.append(self.get(id=uid) )
            except:
                pass
        return allres

    def get( self, id=None, name=None, atID=None ):
        """
        Returns the density of states
        """
        if ( name is None and id is None and atID is None ):
            raise ValueError( "No selection criteria specified!" )

        selectors = [id,name,atID]
        string_selectors = ["id","name","atID"]
        selector_types = {
            "id":"integer",
            "atID":"integer",
            "name":"text"
        }
        for selector,str_selector in zip(selectors,string_selectors):
            if ( selector is None ):
                continue

            conn = sq.connect( self.db_name )
            cur = conn.cursor()
            if ( selector_types[str_selector] == "text" ):
                selector = "\'"+selector+"\'"
            sql = "select id,name,atID,omega_e,dos_e from phonon_dos where {}={}".format(str_selector,selector)
            cur.execute(sql)
            res = cur.fetchone()
            conn.close()

            if ( res is None ):
                raise IOError( "Could not find any match for sql: {}".format(sql))
            res_dict = {}

            res_dict["id"] = int(res[0])
            res_dict["name"] = res[1]
            res_dict["atID"] = res[2]
            res_dict["omega_e"] = self.blob_to_array(res[3])
            res_dict["dos_e"] = self.blob_to_array(res[4])
            return res_dict

    def array_to_blob( self, array ):
        """
        Convert an array to a blob
        """
        try:
            data = np.array( array )
        except Exception as exc:
            print (str(exc))
            raise Exception( "Could not convert passed object to a numpy array!" )

        out = io.BytesIO()
        np.save(out,array)
        out.seek(0)
        return sq.Binary(out.read())

    def blob_to_array( self, blob ):
        """
        Converts a blob to a numpy array
        """
        out = io.BytesIO(blob)
        out.seek(0)
        return np.load(out)

    def get_all_atIDs( self ):
        """
        Returns a list with all atom IDs
        """
        all_entries = self.get_all()
        all_atom_ids = [res["atID"] for res in all_entries]
        return all_atom_ids

    def sync_with_ase_db( self, ase_db_name, selection=[("converged","=",1)] ):
        """
        Syncronizes this database with another ASE database.
        Writes all AtomRow objects into this database
        """
        raise NotImplementedError( "This function does not work yet!" )
        all_atom_ids = self.get_all_atIDs()
        ase_db = connect( ase_db_name )
        ase_db_connection_self = connect( self.db_name )
        num_new_structures = 0
        id_converter = {}
        for row in ase_db.select(selection):
            if ( row.id in all_atom_ids ):
                continue
            kvp = row.key_value_pairs
            atoms = ase_db.get_atoms(id=row.id)
            newID = ase_db_connection_self.write( atoms, key_value_pairs=kvp )
            id_converter[newID] = row.id
            # TODO: Write the new atom objects here
            num_new_structures += 1

        # Have to convert the IDs ASE created back to the original ones
        conn = sq.connect( self.db_name )
        cur = conn.cursor()
        sql = "from systems select id"
        cur.execute( sql )
        entries = cur.fetchall()
        for entry in entries:
            if ( entry[0] in id_converter.keys() ):
                sql = "update systems set id=? where id=?"
                cur.execute( sql, (id_converter[entry[0]], entry[0]) )
        conn.commit()
        cur.close()
        conn.close()
        print ( "Added %d new structures to the database"%(num_new_structures) )
