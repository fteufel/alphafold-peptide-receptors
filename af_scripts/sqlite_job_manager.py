'''
A simple sqlite database so multiple worker nodes can track prediction jobs.
'''

import sqlite3
from datetime import datetime

class SQLiteJobManager():
    def __init__(self, db_path: str):
        
        self.db_path = db_path
        # self
        # self.cursor = self.con.cursor()


        con = sqlite3.connect(self.db_path)
        cur = con.cursor()
        with con:
            cmd = 'CREATE TABLE IF NOT EXISTS af_jobs(peptide_name, receptor_name, peptide_seq, receptor_seq, dir, start_date, end_date, status)'
            cur.execute(cmd)
        
        con.close()

    # inputs.append((rec_id, pep_id, rec_seq, pep_seq))
    def try_add_job(self, pep_id: str, rec_id: str, pep_seq: str, rec_seq: str, out_dir: str):
        '''Add a job to the database. Will return False if the job already exists.'''

        con = sqlite3.connect(self.db_path)
        con.isolation_level = 'EXCLUSIVE'
        with con:
            cur = con.cursor()

            cmd = f"SELECT * FROM af_jobs WHERE (peptide_name='{pep_id}') AND (receptor_name='{rec_id}') AND (dir='{out_dir}')"
            res = cur.execute(cmd)

            if len(res.fetchall())>0:
                # con.close()
                return False
            else:
                cmd = f"INSERT INTO af_jobs VALUES ('{pep_id}', '{rec_id}', '{pep_seq}', '{rec_seq}', '{out_dir}', 'NaN', 'NaN', 'WAITING')"
                cur.execute(cmd)
                con.commit()
                # con.close()

                return True




    def get_job(self):
        '''Retrieve a WAITING job from the database, return it and mark it as RUNNING'''
        con = sqlite3.connect(self.db_path)
        con.isolation_level = 'EXCLUSIVE'
        with con:
            cur = con.cursor()
        
            cmd = f"SELECT * FROM af_jobs WHERE status='WAITING'"
            res = cur.execute(cmd)
            pep_id, rec_id, pep_seq, rec_seq, out_dir, start_date, end_date, status = res.fetchone()

            time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            cmd = f"UPDATE af_jobs SET start_date = '{time}', status='RUNNING' WHERE (peptide_name='{pep_id}') AND (receptor_name='{rec_id}')"
            cur.execute(cmd)
            con.commit()

        con.close()

        return pep_id, rec_id, pep_seq, rec_seq, out_dir


    def mark_job_as_complete(self, pep_id, rec_id, out_dir):
        '''Mark the specified job as COMPLETE'''
        con = sqlite3.connect(self.db_path)
        cur = con.cursor()
        
        time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        cmd = f"UPDATE af_jobs SET end_date = '{time}', status='COMPLETE' WHERE (peptide_name='{pep_id}') AND (receptor_name='{rec_id}') AND (dir='{out_dir}')"
        cur.execute(cmd)
        con.commit()
        con.close()


    def mark_job_as_failed(self, pep_id, rec_id, out_dir):
        '''Mark the specified job as COMPLETE'''
        con = sqlite3.connect(self.db_path)
        cur = con.cursor()
        
        time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        cmd = f"UPDATE af_jobs SET end_date = '{time}', status='FAILED' WHERE (peptide_name='{pep_id}') AND (receptor_name='{rec_id}') AND (dir='{out_dir}')"
        cur.execute(cmd)
        con.commit()
        con.close()

    def reset_failed_jobs(self):
        '''Reset all FAILED jobs to WAITING'''
        con = sqlite3.connect(self.db_path)
        cur = con.cursor()

        cmd = f"UPDATE af_jobs SET start_data = 'NaN', end_date = 'NaN', status='WAITING' WHERE (status='FAILED')"
        cur.execute(cmd)
        con.commit()
        con.close()
