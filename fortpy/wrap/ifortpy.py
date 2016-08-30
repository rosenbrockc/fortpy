import _ifortpy
import f90wrap.runtime
import logging

class Fortpy(f90wrap.runtime.FortranModule):
    """
    Module fortpy
    
    
    Defined at fortpy.f90 lines 4-16128
    
    """
    @staticmethod
    def fpy_set_verbosity(v):
        """
        fpy_set_verbosity(v)
        
        
        Defined at fortpy.f90 lines 163-167
        
        Parameters
        ----------
        v : int
        
        """
        _ifortpy.f90wrap_fpy_set_verbosity(v=v)
    
    @staticmethod
    def fpy_period_join_indices(indices, n):
        """
        pslist = fpy_period_join_indices(indices, n)
        
        
        Defined at fortpy.f90 lines 14214-14237
        
        Parameters
        ----------
        indices : int array
        n : int
        
        Returns
        -------
        pslist : str
        
        """
        pslist = _ifortpy.f90wrap_fpy_period_join_indices(indices=indices, n=n)
        return pslist
    
    @staticmethod
    def random_real(min_bn, max_bn):
        """
        variable = random_real(min_bn, max_bn)
        
        
        Defined at fortpy.f90 lines 14239-14247
        
        Parameters
        ----------
        min_bn : float
        max_bn : float
        
        Returns
        -------
        variable : float
        
        """
        variable = _ifortpy.f90wrap_random_real(min_bn=min_bn, max_bn=max_bn)
        return variable
    
    @staticmethod
    def random_integer(min_bn, max_bn):
        """
        variable = random_integer(min_bn, max_bn)
        
        
        Defined at fortpy.f90 lines 14249-14257
        
        Parameters
        ----------
        min_bn : int
        max_bn : int
        
        Returns
        -------
        variable : int
        
        """
        variable = _ifortpy.f90wrap_random_integer(min_bn=min_bn, max_bn=max_bn)
        return variable
    
    @staticmethod
    def random_init(seed=None):
        """
        random_init([seed])
        
        
        Defined at fortpy.f90 lines 14260-14270
        
        Parameters
        ----------
        seed : int array
        
        """
        _ifortpy.f90wrap_random_init(seed=seed)
    
    @staticmethod
    def char_escape_word(word, escaped):
        """
        char_escape_word(word, escaped)
        
        
        Defined at fortpy.f90 lines 14274-14291
        
        Parameters
        ----------
        word : str
        escaped : str
        
        """
        _ifortpy.f90wrap_char_escape_word(word=word, escaped=escaped)
    
    @staticmethod
    def char_write_trimmed(variable):
        """
        char_write_trimmed(variable)
        
        
        Defined at fortpy.f90 lines 14295-14321
        
        Parameters
        ----------
        variable : str array
        
        """
        _ifortpy.f90wrap_char_write_trimmed(variable=variable)
    
    @staticmethod
    def file_open(filename, n, template_name):
        """
        file_open(filename, n, template_name)
        
        
        Defined at fortpy.f90 lines 15841-15867
        
        Parameters
        ----------
        filename : str
        n : int
        template_name : str
        
        """
        _ifortpy.f90wrap_file_open(filename=filename, n=n, template_name=template_name)
    
    @staticmethod
    def file_close():
        """
        file_close()
        
        
        Defined at fortpy.f90 lines 15869-15873
        
        
        """
        _ifortpy.f90wrap_file_close()
    
    @staticmethod
    def fpy_value_count(line, length, ischar=None):
        """
        fpy_value_count = fpy_value_count(line, length[, ischar])
        
        
        Defined at fortpy.f90 lines 15881-15998
        
        Parameters
        ----------
        line : str
        length : int
        ischar : bool
        
        Returns
        -------
        fpy_value_count : int
        
        """
        fpy_value_count = _ifortpy.f90wrap_fpy_value_count(line=line, length=length, \
            ischar=ischar)
        return fpy_value_count
    
    @staticmethod
    def fpy_linevalue_count(filename, commentchar, ischar=None):
        """
        nlines, nvalues = fpy_linevalue_count(filename, commentchar[, ischar])
        
        
        Defined at fortpy.f90 lines 16053-16109
        
        Parameters
        ----------
        filename : str
        commentchar : str
        ischar : bool
        
        Returns
        -------
        nlines : int
        nvalues : int
        
        """
        nlines, nvalues = _ifortpy.f90wrap_fpy_linevalue_count(filename=filename, \
            commentchar=commentchar, ischar=ischar)
        return nlines, nvalues
    
    @staticmethod
    def fpy_newunit():
        """
        fpy_newunit, unit = fpy_newunit()
        
        
        Defined at fortpy.f90 lines 16113-16128
        
        
        Returns
        -------
        fpy_newunit : int
        unit : int
        
        """
        fpy_newunit, unit = _ifortpy.f90wrap_fpy_newunit()
        return fpy_newunit, unit
    
    @staticmethod
    def _pysave_realsp_0d(variable, filename):
        """
        _pysave_realsp_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14323-14333
        
        Parameters
        ----------
        variable : float
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realsp_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realsp_1d(variable, filename):
        """
        _pysave_realsp_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14335-14348
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realsp_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realsp_2d(variable, filename):
        """
        _pysave_realsp_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14350-14365
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realsp_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realsp_3d(variable, filename):
        """
        _pysave_realsp_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14367-14386
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realsp_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realsp_4d(variable, filename):
        """
        _pysave_realsp_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14388-14409
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realsp_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realsp_5d(variable, filename):
        """
        _pysave_realsp_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14411-14434
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realsp_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realsp_6d(variable, filename):
        """
        _pysave_realsp_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14436-14461
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realsp_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realsp_7d(variable, filename):
        """
        _pysave_realsp_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14463-14489
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realsp_7d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realdp_0d(variable, filename):
        """
        _pysave_realdp_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14491-14501
        
        Parameters
        ----------
        variable : float
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realdp_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realdp_1d(variable, filename):
        """
        _pysave_realdp_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14503-14516
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realdp_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realdp_2d(variable, filename):
        """
        _pysave_realdp_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14518-14533
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realdp_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realdp_3d(variable, filename):
        """
        _pysave_realdp_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14535-14554
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realdp_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realdp_4d(variable, filename):
        """
        _pysave_realdp_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14556-14577
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realdp_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realdp_5d(variable, filename):
        """
        _pysave_realdp_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14579-14602
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realdp_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realdp_6d(variable, filename):
        """
        _pysave_realdp_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14604-14629
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realdp_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_realdp_7d(variable, filename):
        """
        _pysave_realdp_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14631-14658
        
        Parameters
        ----------
        variable : float array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_realdp_7d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integer_0d(variable, filename):
        """
        _pysave_integer_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14660-14670
        
        Parameters
        ----------
        variable : int
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integer_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integer_1d(variable, filename):
        """
        _pysave_integer_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14672-14685
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integer_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integer_2d(variable, filename):
        """
        _pysave_integer_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14687-14702
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integer_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integer_3d(variable, filename):
        """
        _pysave_integer_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14704-14723
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integer_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integer_4d(variable, filename):
        """
        _pysave_integer_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14725-14746
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integer_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integer_5d(variable, filename):
        """
        _pysave_integer_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14748-14771
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integer_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integer_6d(variable, filename):
        """
        _pysave_integer_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14773-14798
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integer_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integer_7d(variable, filename):
        """
        _pysave_integer_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14800-14826
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integer_7d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integersp_0d(variable, filename):
        """
        _pysave_integersp_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14828-14838
        
        Parameters
        ----------
        variable : int
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integersp_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integersp_1d(variable, filename):
        """
        _pysave_integersp_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14840-14853
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integersp_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integersp_2d(variable, filename):
        """
        _pysave_integersp_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14855-14870
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integersp_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integersp_3d(variable, filename):
        """
        _pysave_integersp_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14872-14891
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integersp_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integersp_4d(variable, filename):
        """
        _pysave_integersp_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14893-14914
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integersp_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integersp_5d(variable, filename):
        """
        _pysave_integersp_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14916-14939
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integersp_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integersp_6d(variable, filename):
        """
        _pysave_integersp_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14941-14966
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integersp_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integersp_7d(variable, filename):
        """
        _pysave_integersp_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14968-14994
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integersp_7d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integerdp_0d(variable, filename):
        """
        _pysave_integerdp_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 14996-15006
        
        Parameters
        ----------
        variable : int
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integerdp_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integerdp_1d(variable, filename):
        """
        _pysave_integerdp_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15008-15021
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integerdp_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integerdp_2d(variable, filename):
        """
        _pysave_integerdp_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15023-15038
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integerdp_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integerdp_3d(variable, filename):
        """
        _pysave_integerdp_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15040-15059
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integerdp_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integerdp_4d(variable, filename):
        """
        _pysave_integerdp_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15061-15082
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integerdp_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integerdp_5d(variable, filename):
        """
        _pysave_integerdp_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15084-15107
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integerdp_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integerdp_6d(variable, filename):
        """
        _pysave_integerdp_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15109-15134
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integerdp_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_integerdp_7d(variable, filename):
        """
        _pysave_integerdp_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15136-15163
        
        Parameters
        ----------
        variable : int array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_integerdp_7d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexsp_0d(variable, filename):
        """
        _pysave_complexsp_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15165-15175
        
        Parameters
        ----------
        variable : complex
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexsp_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexsp_1d(variable, filename):
        """
        _pysave_complexsp_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15177-15190
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexsp_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexsp_2d(variable, filename):
        """
        _pysave_complexsp_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15192-15207
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexsp_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexsp_3d(variable, filename):
        """
        _pysave_complexsp_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15209-15228
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexsp_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexsp_4d(variable, filename):
        """
        _pysave_complexsp_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15230-15251
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexsp_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexsp_5d(variable, filename):
        """
        _pysave_complexsp_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15253-15276
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexsp_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexsp_6d(variable, filename):
        """
        _pysave_complexsp_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15278-15303
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexsp_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexsp_7d(variable, filename):
        """
        _pysave_complexsp_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15305-15331
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexsp_7d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexdp_0d(variable, filename):
        """
        _pysave_complexdp_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15333-15343
        
        Parameters
        ----------
        variable : complex
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexdp_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexdp_1d(variable, filename):
        """
        _pysave_complexdp_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15345-15358
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexdp_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexdp_2d(variable, filename):
        """
        _pysave_complexdp_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15360-15375
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexdp_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexdp_3d(variable, filename):
        """
        _pysave_complexdp_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15377-15396
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexdp_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexdp_4d(variable, filename):
        """
        _pysave_complexdp_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15398-15419
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexdp_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexdp_5d(variable, filename):
        """
        _pysave_complexdp_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15421-15444
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexdp_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexdp_6d(variable, filename):
        """
        _pysave_complexdp_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15446-15471
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexdp_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_complexdp_7d(variable, filename):
        """
        _pysave_complexdp_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15473-15500
        
        Parameters
        ----------
        variable : complex array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_complexdp_7d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_character_0d(variable, filename):
        """
        _pysave_character_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15502-15512
        
        Parameters
        ----------
        variable : str
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_character_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_character_1d(variable, filename):
        """
        _pysave_character_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15514-15527
        
        Parameters
        ----------
        variable : str array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_character_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_character_2d(variable, filename):
        """
        _pysave_character_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15529-15544
        
        Parameters
        ----------
        variable : str array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_character_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_character_3d(variable, filename):
        """
        _pysave_character_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15546-15565
        
        Parameters
        ----------
        variable : str array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_character_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_character_4d(variable, filename):
        """
        _pysave_character_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15567-15588
        
        Parameters
        ----------
        variable : str array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_character_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_character_5d(variable, filename):
        """
        _pysave_character_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15590-15613
        
        Parameters
        ----------
        variable : str array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_character_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_character_6d(variable, filename):
        """
        _pysave_character_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15615-15640
        
        Parameters
        ----------
        variable : str array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_character_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_character_7d(variable, filename):
        """
        _pysave_character_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15642-15669
        
        Parameters
        ----------
        variable : str array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_character_7d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_logical_0d(variable, filename):
        """
        _pysave_logical_0d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15671-15681
        
        Parameters
        ----------
        variable : bool
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_logical_0d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_logical_1d(variable, filename):
        """
        _pysave_logical_1d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15683-15696
        
        Parameters
        ----------
        variable : bool array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_logical_1d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_logical_2d(variable, filename):
        """
        _pysave_logical_2d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15698-15713
        
        Parameters
        ----------
        variable : bool array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_logical_2d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_logical_3d(variable, filename):
        """
        _pysave_logical_3d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15715-15734
        
        Parameters
        ----------
        variable : bool array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_logical_3d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_logical_4d(variable, filename):
        """
        _pysave_logical_4d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15736-15757
        
        Parameters
        ----------
        variable : bool array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_logical_4d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_logical_5d(variable, filename):
        """
        _pysave_logical_5d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15759-15782
        
        Parameters
        ----------
        variable : bool array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_logical_5d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_logical_6d(variable, filename):
        """
        _pysave_logical_6d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15784-15809
        
        Parameters
        ----------
        variable : bool array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_logical_6d(variable=variable, filename=filename)
    
    @staticmethod
    def _pysave_logical_7d(variable, filename):
        """
        _pysave_logical_7d(variable, filename)
        
        
        Defined at fortpy.f90 lines 15811-15839
        
        Parameters
        ----------
        variable : bool array
        filename : str
        
        """
        _ifortpy.f90wrap_pysave_logical_7d(variable=variable, filename=filename)
    
    @staticmethod
    def pysave(*args, **kwargs):
        """
        pysave(*args, **kwargs)
        
        
        Defined at fortpy.f90 lines 27-52
        
        Overloaded interface containing the following procedures:
          _pysave_realsp_0d
          _pysave_realsp_1d
          _pysave_realsp_2d
          _pysave_realsp_3d
          _pysave_realsp_4d
          _pysave_realsp_5d
          _pysave_realsp_6d
          _pysave_realsp_7d
          _pysave_realdp_0d
          _pysave_realdp_1d
          _pysave_realdp_2d
          _pysave_realdp_3d
          _pysave_realdp_4d
          _pysave_realdp_5d
          _pysave_realdp_6d
          _pysave_realdp_7d
          _pysave_integer_0d
          _pysave_integer_1d
          _pysave_integer_2d
          _pysave_integer_3d
          _pysave_integer_4d
          _pysave_integer_5d
          _pysave_integer_6d
          _pysave_integer_7d
          _pysave_integersp_0d
          _pysave_integersp_1d
          _pysave_integersp_2d
          _pysave_integersp_3d
          _pysave_integersp_4d
          _pysave_integersp_5d
          _pysave_integersp_6d
          _pysave_integersp_7d
          _pysave_integerdp_0d
          _pysave_integerdp_1d
          _pysave_integerdp_2d
          _pysave_integerdp_3d
          _pysave_integerdp_4d
          _pysave_integerdp_5d
          _pysave_integerdp_6d
          _pysave_integerdp_7d
          _pysave_complexsp_0d
          _pysave_complexsp_1d
          _pysave_complexsp_2d
          _pysave_complexsp_3d
          _pysave_complexsp_4d
          _pysave_complexsp_5d
          _pysave_complexsp_6d
          _pysave_complexsp_7d
          _pysave_complexdp_0d
          _pysave_complexdp_1d
          _pysave_complexdp_2d
          _pysave_complexdp_3d
          _pysave_complexdp_4d
          _pysave_complexdp_5d
          _pysave_complexdp_6d
          _pysave_complexdp_7d
          _pysave_character_0d
          _pysave_character_1d
          _pysave_character_2d
          _pysave_character_3d
          _pysave_character_4d
          _pysave_character_5d
          _pysave_character_6d
          _pysave_character_7d
          _pysave_logical_0d
          _pysave_logical_1d
          _pysave_logical_2d
          _pysave_logical_3d
          _pysave_logical_4d
          _pysave_logical_5d
          _pysave_logical_6d
          _pysave_logical_7d
        
        """
        for proc in [Fortpy._pysave_realsp_0d, Fortpy._pysave_realsp_1d, \
            Fortpy._pysave_realsp_2d, Fortpy._pysave_realsp_3d, \
            Fortpy._pysave_realsp_4d, Fortpy._pysave_realsp_5d, \
            Fortpy._pysave_realsp_6d, Fortpy._pysave_realsp_7d, \
            Fortpy._pysave_realdp_0d, Fortpy._pysave_realdp_1d, \
            Fortpy._pysave_realdp_2d, Fortpy._pysave_realdp_3d, \
            Fortpy._pysave_realdp_4d, Fortpy._pysave_realdp_5d, \
            Fortpy._pysave_realdp_6d, Fortpy._pysave_realdp_7d, \
            Fortpy._pysave_integer_0d, Fortpy._pysave_integer_1d, \
            Fortpy._pysave_integer_2d, Fortpy._pysave_integer_3d, \
            Fortpy._pysave_integer_4d, Fortpy._pysave_integer_5d, \
            Fortpy._pysave_integer_6d, Fortpy._pysave_integer_7d, \
            Fortpy._pysave_integersp_0d, Fortpy._pysave_integersp_1d, \
            Fortpy._pysave_integersp_2d, Fortpy._pysave_integersp_3d, \
            Fortpy._pysave_integersp_4d, Fortpy._pysave_integersp_5d, \
            Fortpy._pysave_integersp_6d, Fortpy._pysave_integersp_7d, \
            Fortpy._pysave_integerdp_0d, Fortpy._pysave_integerdp_1d, \
            Fortpy._pysave_integerdp_2d, Fortpy._pysave_integerdp_3d, \
            Fortpy._pysave_integerdp_4d, Fortpy._pysave_integerdp_5d, \
            Fortpy._pysave_integerdp_6d, Fortpy._pysave_integerdp_7d, \
            Fortpy._pysave_complexsp_0d, Fortpy._pysave_complexsp_1d, \
            Fortpy._pysave_complexsp_2d, Fortpy._pysave_complexsp_3d, \
            Fortpy._pysave_complexsp_4d, Fortpy._pysave_complexsp_5d, \
            Fortpy._pysave_complexsp_6d, Fortpy._pysave_complexsp_7d, \
            Fortpy._pysave_complexdp_0d, Fortpy._pysave_complexdp_1d, \
            Fortpy._pysave_complexdp_2d, Fortpy._pysave_complexdp_3d, \
            Fortpy._pysave_complexdp_4d, Fortpy._pysave_complexdp_5d, \
            Fortpy._pysave_complexdp_6d, Fortpy._pysave_complexdp_7d, \
            Fortpy._pysave_character_0d, Fortpy._pysave_character_1d, \
            Fortpy._pysave_character_2d, Fortpy._pysave_character_3d, \
            Fortpy._pysave_character_4d, Fortpy._pysave_character_5d, \
            Fortpy._pysave_character_6d, Fortpy._pysave_character_7d, \
            Fortpy._pysave_logical_0d, Fortpy._pysave_logical_1d, \
            Fortpy._pysave_logical_2d, Fortpy._pysave_logical_3d, \
            Fortpy._pysave_logical_4d, Fortpy._pysave_logical_5d, \
            Fortpy._pysave_logical_6d, Fortpy._pysave_logical_7d]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _fpy_read_realsp_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realsp_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 222-276
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realsp_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realdp_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realdp_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 774-828
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realdp_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integer_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integer_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 1327-1381
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integer_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integersp_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integersp_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 1879-1933
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integersp_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integerdp_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integerdp_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 2431-2485
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integerdp_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexsp_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexsp_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 2984-3038
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexsp_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexdp_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexdp_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 3536-3590
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexdp_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_character_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_character_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 4089-4143
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : str
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_character_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_logical_0d(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_logical_0d(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 4642-4696
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : bool
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_logical_0d(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def fpy_read(*args, **kwargs):
        """
        fpy_read(*args, **kwargs)
        
        
        Defined at fortpy.f90 lines 55-80
        
        Overloaded interface containing the following procedures:
          _fpy_read_realsp_0d
          _fpy_read_realdp_0d
          _fpy_read_integer_0d
          _fpy_read_integersp_0d
          _fpy_read_integerdp_0d
          _fpy_read_complexsp_0d
          _fpy_read_complexdp_0d
          _fpy_read_character_0d
          _fpy_read_logical_0d
        
        """
        for proc in [Fortpy._fpy_read_realsp_0d, Fortpy._fpy_read_realdp_0d, \
            Fortpy._fpy_read_integer_0d, Fortpy._fpy_read_integersp_0d, \
            Fortpy._fpy_read_integerdp_0d, Fortpy._fpy_read_complexsp_0d, \
            Fortpy._fpy_read_complexdp_0d, Fortpy._fpy_read_character_0d, \
            Fortpy._fpy_read_logical_0d]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @staticmethod
    def _fpy_read_realsp_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realsp_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5196-5264
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realsp_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realsp_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realsp_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5266-5322
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realsp_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realsp_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realsp_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5324-5394
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realsp_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realsp_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realsp_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5396-5467
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realsp_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realsp_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realsp_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5469-5541
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realsp_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realsp_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realsp_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5543-5616
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realsp_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realsp_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realsp_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5618-5691
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realsp_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realdp_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realdp_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5693-5761
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realdp_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realdp_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realdp_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5763-5819
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realdp_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realdp_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realdp_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5821-5891
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realdp_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realdp_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realdp_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5893-5964
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realdp_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realdp_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realdp_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 5966-6038
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realdp_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realdp_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realdp_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6040-6113
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realdp_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_realdp_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_realdp_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6115-6189
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : float array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_realdp_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integer_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integer_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6191-6259
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integer_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integer_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integer_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6261-6317
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integer_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integer_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integer_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6319-6389
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integer_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integer_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integer_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6391-6462
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integer_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integer_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integer_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6464-6536
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integer_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integer_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integer_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6538-6611
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integer_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integer_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integer_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6613-6686
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integer_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integersp_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integersp_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6688-6756
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integersp_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integersp_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integersp_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6758-6814
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integersp_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integersp_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integersp_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6816-6886
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integersp_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integersp_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integersp_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6888-6959
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integersp_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integersp_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integersp_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 6961-7033
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integersp_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integersp_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integersp_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7035-7108
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integersp_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integersp_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integersp_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7110-7183
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integersp_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integerdp_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integerdp_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7185-7253
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integerdp_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integerdp_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integerdp_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7255-7311
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integerdp_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integerdp_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integerdp_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7313-7383
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integerdp_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integerdp_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integerdp_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7385-7456
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integerdp_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integerdp_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integerdp_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7458-7530
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integerdp_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integerdp_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integerdp_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7532-7605
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integerdp_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_integerdp_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_integerdp_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7607-7681
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : int array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_integerdp_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexsp_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexsp_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7683-7751
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexsp_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexsp_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexsp_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7753-7809
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexsp_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexsp_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexsp_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7811-7881
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexsp_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexsp_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexsp_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7883-7954
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexsp_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexsp_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexsp_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 7956-8028
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexsp_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexsp_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexsp_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8030-8103
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexsp_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexsp_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexsp_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8105-8178
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexsp_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexdp_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexdp_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8180-8248
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexdp_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexdp_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexdp_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8250-8306
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexdp_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexdp_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexdp_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8308-8378
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexdp_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexdp_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexdp_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8380-8451
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexdp_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexdp_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexdp_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8453-8525
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexdp_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexdp_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexdp_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8527-8600
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexdp_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_complexdp_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_complexdp_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8602-8676
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : complex array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_complexdp_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_character_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_character_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8678-8746
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : str array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_character_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_character_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_character_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8748-8804
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : str array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_character_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_character_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_character_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8806-8876
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : str array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_character_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_character_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_character_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8878-8949
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : str array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_character_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_character_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_character_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 8951-9023
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : str array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_character_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_character_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_character_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9025-9098
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : str array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_character_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_character_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_character_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9100-9174
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : str array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_character_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_logical_1df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_logical_1df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9176-9244
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : bool array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_logical_1df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_logical_2df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_logical_2df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9246-9302
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : bool array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_logical_2df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_logical_3df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_logical_3df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9304-9374
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : bool array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_logical_3df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_logical_4df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_logical_4df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9376-9447
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : bool array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_logical_4df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_logical_5df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_logical_5df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9449-9521
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : bool array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_logical_5df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_logical_6df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_logical_6df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9523-9596
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : bool array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_logical_6df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def _fpy_read_logical_7df(filename, commentchar, variable, strict_=None):
        """
        success_ = _fpy_read_logical_7df(filename, commentchar, variable[, strict_])
        
        
        Defined at fortpy.f90 lines 9598-9673
        
        Parameters
        ----------
        filename : str
        commentchar : str
        variable : bool array
        strict_ : bool
        
        Returns
        -------
        success_ : bool
        
        """
        success_ = _ifortpy.f90wrap_fpy_read_logical_7df(filename=filename, \
            commentchar=commentchar, variable=variable, strict_=strict_)
        return success_
    
    @staticmethod
    def fpy_read_f(*args, **kwargs):
        """
        fpy_read_f(*args, **kwargs)
        
        
        Defined at fortpy.f90 lines 84-112
        
        Overloaded interface containing the following procedures:
          _fpy_read_realsp_1df
          _fpy_read_realsp_2df
          _fpy_read_realsp_3df
          _fpy_read_realsp_4df
          _fpy_read_realsp_5df
          _fpy_read_realsp_6df
          _fpy_read_realsp_7df
          _fpy_read_realdp_1df
          _fpy_read_realdp_2df
          _fpy_read_realdp_3df
          _fpy_read_realdp_4df
          _fpy_read_realdp_5df
          _fpy_read_realdp_6df
          _fpy_read_realdp_7df
          _fpy_read_integer_1df
          _fpy_read_integer_2df
          _fpy_read_integer_3df
          _fpy_read_integer_4df
          _fpy_read_integer_5df
          _fpy_read_integer_6df
          _fpy_read_integer_7df
          _fpy_read_integersp_1df
          _fpy_read_integersp_2df
          _fpy_read_integersp_3df
          _fpy_read_integersp_4df
          _fpy_read_integersp_5df
          _fpy_read_integersp_6df
          _fpy_read_integersp_7df
          _fpy_read_integerdp_1df
          _fpy_read_integerdp_2df
          _fpy_read_integerdp_3df
          _fpy_read_integerdp_4df
          _fpy_read_integerdp_5df
          _fpy_read_integerdp_6df
          _fpy_read_integerdp_7df
          _fpy_read_complexsp_1df
          _fpy_read_complexsp_2df
          _fpy_read_complexsp_3df
          _fpy_read_complexsp_4df
          _fpy_read_complexsp_5df
          _fpy_read_complexsp_6df
          _fpy_read_complexsp_7df
          _fpy_read_complexdp_1df
          _fpy_read_complexdp_2df
          _fpy_read_complexdp_3df
          _fpy_read_complexdp_4df
          _fpy_read_complexdp_5df
          _fpy_read_complexdp_6df
          _fpy_read_complexdp_7df
          _fpy_read_character_1df
          _fpy_read_character_2df
          _fpy_read_character_3df
          _fpy_read_character_4df
          _fpy_read_character_5df
          _fpy_read_character_6df
          _fpy_read_character_7df
          _fpy_read_logical_1df
          _fpy_read_logical_2df
          _fpy_read_logical_3df
          _fpy_read_logical_4df
          _fpy_read_logical_5df
          _fpy_read_logical_6df
          _fpy_read_logical_7df
        
        """
        for proc in [Fortpy._fpy_read_realsp_1df, Fortpy._fpy_read_realsp_2df, \
            Fortpy._fpy_read_realsp_3df, Fortpy._fpy_read_realsp_4df, \
            Fortpy._fpy_read_realsp_5df, Fortpy._fpy_read_realsp_6df, \
            Fortpy._fpy_read_realsp_7df, Fortpy._fpy_read_realdp_1df, \
            Fortpy._fpy_read_realdp_2df, Fortpy._fpy_read_realdp_3df, \
            Fortpy._fpy_read_realdp_4df, Fortpy._fpy_read_realdp_5df, \
            Fortpy._fpy_read_realdp_6df, Fortpy._fpy_read_realdp_7df, \
            Fortpy._fpy_read_integer_1df, Fortpy._fpy_read_integer_2df, \
            Fortpy._fpy_read_integer_3df, Fortpy._fpy_read_integer_4df, \
            Fortpy._fpy_read_integer_5df, Fortpy._fpy_read_integer_6df, \
            Fortpy._fpy_read_integer_7df, Fortpy._fpy_read_integersp_1df, \
            Fortpy._fpy_read_integersp_2df, Fortpy._fpy_read_integersp_3df, \
            Fortpy._fpy_read_integersp_4df, Fortpy._fpy_read_integersp_5df, \
            Fortpy._fpy_read_integersp_6df, Fortpy._fpy_read_integersp_7df, \
            Fortpy._fpy_read_integerdp_1df, Fortpy._fpy_read_integerdp_2df, \
            Fortpy._fpy_read_integerdp_3df, Fortpy._fpy_read_integerdp_4df, \
            Fortpy._fpy_read_integerdp_5df, Fortpy._fpy_read_integerdp_6df, \
            Fortpy._fpy_read_integerdp_7df, Fortpy._fpy_read_complexsp_1df, \
            Fortpy._fpy_read_complexsp_2df, Fortpy._fpy_read_complexsp_3df, \
            Fortpy._fpy_read_complexsp_4df, Fortpy._fpy_read_complexsp_5df, \
            Fortpy._fpy_read_complexsp_6df, Fortpy._fpy_read_complexsp_7df, \
            Fortpy._fpy_read_complexdp_1df, Fortpy._fpy_read_complexdp_2df, \
            Fortpy._fpy_read_complexdp_3df, Fortpy._fpy_read_complexdp_4df, \
            Fortpy._fpy_read_complexdp_5df, Fortpy._fpy_read_complexdp_6df, \
            Fortpy._fpy_read_complexdp_7df, Fortpy._fpy_read_character_1df, \
            Fortpy._fpy_read_character_2df, Fortpy._fpy_read_character_3df, \
            Fortpy._fpy_read_character_4df, Fortpy._fpy_read_character_5df, \
            Fortpy._fpy_read_character_6df, Fortpy._fpy_read_character_7df, \
            Fortpy._fpy_read_logical_1df, Fortpy._fpy_read_logical_2df, \
            Fortpy._fpy_read_logical_3df, Fortpy._fpy_read_logical_4df, \
            Fortpy._fpy_read_logical_5df, Fortpy._fpy_read_logical_6df, \
            Fortpy._fpy_read_logical_7df]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    @property
    def fileunit(self):
        """
        Element fileunit ftype=integer  pytype=int
        
        
        Defined at fortpy.f90 line 11
        
        """
        return _ifortpy.f90wrap_fortpy__get__fileunit()
    
    @fileunit.setter
    def fileunit(self, fileunit):
        _ifortpy.f90wrap_fortpy__set__fileunit(fileunit)
    
    @property
    def fpy_verbose(self):
        """
        Element fpy_verbose ftype=integer  pytype=int
        
        
        Defined at fortpy.f90 line 13
        
        """
        return _ifortpy.f90wrap_fortpy__get__fpy_verbose()
    
    @fpy_verbose.setter
    def fpy_verbose(self, fpy_verbose):
        _ifortpy.f90wrap_fortpy__set__fpy_verbose(fpy_verbose)
    
    @property
    def seeded(self):
        """
        Element seeded ftype=logical pytype=bool
        
        
        Defined at fortpy.f90 line 17
        
        """
        return _ifortpy.f90wrap_fortpy__get__seeded()
    
    @seeded.setter
    def seeded(self, seeded):
        _ifortpy.f90wrap_fortpy__set__seeded(seeded)
    
    @property
    def fdp(self):
        """
        Element fdp ftype=integer pytype=int
        
        
        Defined at fortpy.f90 line 19
        
        """
        return _ifortpy.f90wrap_fortpy__get__fdp()
    
    @property
    def fsp(self):
        """
        Element fsp ftype=integer pytype=int
        
        
        Defined at fortpy.f90 line 20
        
        """
        return _ifortpy.f90wrap_fortpy__get__fsp()
    
    @property
    def fsi(self):
        """
        Element fsi ftype=integer pytype=int
        
        
        Defined at fortpy.f90 line 21
        
        """
        return _ifortpy.f90wrap_fortpy__get__fsi()
    
    @property
    def fli(self):
        """
        Element fli ftype=integer pytype=int
        
        
        Defined at fortpy.f90 line 22
        
        """
        return _ifortpy.f90wrap_fortpy__get__fli()
    
    def __str__(self):
        ret = ['<fortpy>{\n']
        ret.append('    fileunit : ')
        ret.append(repr(self.fileunit))
        ret.append(',\n    fpy_verbose : ')
        ret.append(repr(self.fpy_verbose))
        ret.append(',\n    seeded : ')
        ret.append(repr(self.seeded))
        ret.append(',\n    fdp : ')
        ret.append(repr(self.fdp))
        ret.append(',\n    fsp : ')
        ret.append(repr(self.fsp))
        ret.append(',\n    fsi : ')
        ret.append(repr(self.fsi))
        ret.append(',\n    fli : ')
        ret.append(repr(self.fli))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

fortpy = Fortpy()

