n_sources = 3
n_sinks = 3
m = -0.9
seqsource_type = 'NUCL_U_POL'
i = 0
t_sink = 14
s_size = 16
t_size = 64

def smear_factor(x):
    return 1+x

def smear_iter(x):
    return 5+5*x

f = open('bar3ptfnSEED.ini.xml', 'w')

f.write(
        '<?xml version="1.0"?>\n'+
        '<chroma>\n'+
        '<annotation>\n'+
        'Calculate bar3ptfn - insertions in 3-pt functions\n'+
        '</annotation>\n'+
        '<Param> \n'+
        '  <InlineMeasurements>\n\n')

for x in range(n_sources):
    f.write('    <elem>\n'+
            '      <Name>MAKE_SOURCE</Name>\n'+
            '      <Frequency>1</Frequency>\n'+
            '      <Param>\n'+
            '        <version>6</version>\n'+
            '        <Source>\n'+
            '          <version>1</version>\n'+
            '          <SourceType>SHELL_SOURCE</SourceType>\n'+
            '          <j_decay>3</j_decay>\n'+
            '          <t_srce>0 0 0 0</t_srce>\n'+
            '          <SmearingParam>\n'+
            '              <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>\n'+
            '              <wvf_param>{0}</wvf_param>\n'.format(smear_factor(x))+
            '              <wvfIntPar>{0}</wvfIntPar>\n'.format(smear_iter(x))+
            '              <no_smear_dir>3</no_smear_dir>\n'+
            '            </SmearingParam>\n'+
            '      <Displacement>\n'+
            '              <version>1</version>\n'+
            '              <DisplacementType>NONE</DisplacementType>\n'+
            '            </Displacement>\n'+
            '          <LinkSmearing>\n'+
            '              <LinkSmearingType>APE_SMEAR</LinkSmearingType>\n'+
            '              <link_smear_fact>{0}</link_smear_fact>\n'.format(smear_factor(x))+
            '              <link_smear_num>{0}</link_smear_num>\n'.format(smear_iter(x))+
            '              <no_smear_dir>3</no_smear_dir>\n'+
            '            </LinkSmearing>\n'+
            '      </Source>\n'+
            '      </Param>\n'+
            '        <NamedObject>\n'+
            '          <gauge_id>default_gauge_field</gauge_id>\n'+
            '          <source_id>sh_{0}_source</source_id>\n'.format(x)+
            '        </NamedObject>\n'+
            '    </elem>\n\n')
    f.write('   <elem>\n'+
            '      <Name>PROPAGATOR</Name>\n'+
            '      <Frequency>1</Frequency>\n'+
            '      <Param>\n'+
            '        <version>10</version>\n'+
            '        <quarkSpinType>FULL</quarkSpinType>\n'+
            '        <obsvP>false</obsvP>\n'+
            '        <numRetries>1</numRetries>\n'+
            '        <FermionAction>\n'+
            '         <FermAct>WILSON</FermAct>\n'+
            '         <Mass>{0}</Mass>\n'.format(m)+
            '         <AnisoParam>\n'+
            '           <anisoP>true</anisoP>\n'+
            '           <t_dir>3</t_dir>\n'+
            '           <xi_0>2.0</xi_0>\n'+
            '           <nu>2.0</nu>\n'+
            '         </AnisoParam>\n'+
            '         <FermionBC>\n'+
            '           <FermBC>SIMPLE_FERMBC</FermBC>\n'+
            '           <boundary>1 1 1 -1</boundary>\n'+
            '         </FermionBC>\n'+
            '        </FermionAction>\n'+
            '    <InvertParam>\n'+
            '        <invType>QUDA_WILSON_INVERTER</invType>\n'+
            '        <WilsonParams>\n'+
            '            <Mass>{0}</Mass>\n'.format(m)+
            '            <AnisoParam>\n'+
            '              <anisoP>true</anisoP>\n'+
            '              <t_dir>3</t_dir>\n'+
            '              <xi_0>2.0</xi_0>\n'+
            '              <nu>2.0</nu>\n'+
            '            </AnisoParam>\n'+
            '         </WilsonParams>\n'+
            '          <RsdTarget>1e-8</RsdTarget>\n'+
            '          <RsdToleranceFactor>1000</RsdToleranceFactor>\n'+
            '          <Delta>0.5</Delta>\n'+
            '          <MaxIter>10000</MaxIter>\n'+
            '          <AntiPeriodicT>true</AntiPeriodicT>\n'+
            '          <SolverType>CG</SolverType>\n'+
            '          <Verbose>false</Verbose>\n'+
            '          <AsymmetricLinop>false</AsymmetricLinop>\n'+
            '          <CudaPrecision>SINGLE</CudaPrecision>\n'+
            '          <CudaReconstruct>RECONS_NONE</CudaReconstruct>\n'+
            '          <CudaSloppyPrecision>HALF</CudaSloppyPrecision>\n'+
            '          <CudaSloppyReconstruct>RECONS_NONE</CudaSloppyReconstruct>\n'+
            '          <AxialGaugeFix>false</AxialGaugeFix>\n'+
            '          <AutotuneDslash>true</AutotuneDslash>\n'+
            '        </InvertParam>\n'+
            '      </Param>\n'+
            '       <NamedObject>\n'+
            '        <gauge_id>default_gauge_field</gauge_id>\n'+
            '        <source_id>sh_{0}_source</source_id>\n'.format(x)+
            '        <prop_id>sh_{0}_prop</prop_id>\n'.format(x)+
            '       </NamedObject>\n'+
            '    </elem>\n\n')
    f.write('   <elem>\n'+
            '      <Name>ERASE_NAMED_OBJECT</Name>\n'+
            '      <Frequency>1</Frequency>\n'+
            '      <NamedObject>\n'+
            '        <object_id>sh_{0}_source</object_id>\n'.format(x)+
            '      </NamedObject>\n'+
            '    </elem>\n\n')
    for y in range(n_sinks):
        f.write('    <elem>\n'+
                '      <Name>SINK_SMEAR</Name>\n'+
                '      <Frequency>1</Frequency>\n'+
                '      <Param>\n'+
                '        <version>5</version>\n'+
                '        <Sink>\n'+
                '          <version>2</version>\n'+
                '          <SinkType>SHELL_SINK</SinkType>\n'+
                '          <j_decay>3</j_decay>\n'+
                '          <Displacement>\n'+
                '            <version>1</version>\n'+
                '            <DisplacementType>NONE</DisplacementType>\n'+
                '          </Displacement>\n'+
                '        <SmearingParam>\n'+
                '            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>\n'+
                '              <wvf_param>{0}</wvf_param>\n'.format(smear_factor(y))+
                '              <wvfIntPar>{0}</wvfIntPar>\n'.format(smear_iter(y))+
                '            <no_smear_dir>3</no_smear_dir>\n'+
                '          </SmearingParam>\n'+
                '          <LinkSmearing>\n'+
                '            <LinkSmearingType>APE_SMEAR</LinkSmearingType>\n'+
                '            <link_smear_fact>{0}</link_smear_fact>\n'.format(smear_factor(y))+
                '            <link_smear_num>{0}</link_smear_num>\n'.format(smear_iter(y))+
                '            <no_smear_dir>3</no_smear_dir>\n'+
                '          </LinkSmearing>\n'+
                '        </Sink>\n'+
                '      </Param>\n'+
                '      <NamedObject>\n'+
                '        <gauge_id>default_gauge_field</gauge_id>\n'+
                '        <prop_id>sh_{0}_prop</prop_id>\n'.format(x)+
                '        <smeared_prop_id>sh_{0}_sh_{1}_sink</smeared_prop_id>\n'.format(x,y)+
                '      </NamedObject>\n'+
                '    </elem>\n\n')
        f.write('    <elem>\n'+
                '      <Name>HADRON_SPECTRUM</Name>\n'+
                '      <Frequency>1</Frequency>\n'+
                '      <Param>\n'+
                '        <version>1</version>\n'+
                '        <MesonP>true</MesonP>\n'+
                '        <CurrentP>true</CurrentP>\n'+
                '        <BaryonP>true</BaryonP>\n'+
                '        <time_rev>false</time_rev>\n'+
                '        <mom2_max>1</mom2_max>\n'+
                '        <avg_equiv_mom>false</avg_equiv_mom>\n'+
                '       </Param>\n'+
                '      <NamedObject>\n'+
                '        <gauge_id>default_gauge_field</gauge_id>\n'+
                '        <sink_pairs>\n'+
                '          <elem>\n'+
                '            <first_id>sh_{0}_sh_{1}_sink</first_id>\n'.format(x,y)+
                '            <second_id>sh_{0}_sh_{1}_sink</second_id>\n'.format(x,y)+
                '          </elem>\n'+
                '        </sink_pairs>\n'+
                '      </NamedObject>\n'+
                '      <xml_file>hadspec_sh_{0}_sh_{1}_CONFIGNAME.dat.xml</xml_file>\n'.format(x,y)+
                '    </elem>\n\n')
        f.write('   <elem>\n'+
                '      <Name>ERASE_NAMED_OBJECT</Name>\n'+
                '      <Frequency>1</Frequency>\n'+
                '      <NamedObject>\n'+
                '        <object_id>sh_{0}_sh_{1}_sink</object_id>\n'.format(x,y)+
                '      </NamedObject>\n'+
                '    </elem>\n\n')
        f.write(
                '    <elem>\n'+
                '      <Name>SEQSOURCE</Name>\n'+
                '      <Frequency>1</Frequency>\n'+
                '      <Param>\n'+
                '        <version>1</version>\n'+
                '        <seq_src>{0}</seq_src>\n'.format(seqsource_type)+
                '        <t_sink>{0}</t_sink>\n'.format(t_sink)+
                '        <sink_mom>0 0 0</sink_mom>\n'+
                '      </Param>\n'+
                '      <PropSink>\n'+
                '        <version>5</version>\n'+
                '        <Sink>\n'+
                '          <version>2</version>\n'+
                '          <SinkType>SHELL_SINK</SinkType>\n'+
                '          <j_decay>3</j_decay>\n'+
                '          <Displacement>\n'+
                '            <version>1</version>\n'+
                '            <DisplacementType>NONE</DisplacementType>\n'+
                '          </Displacement>\n'+
                '        <SmearingParam>\n'+
                '            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>\n'+
                '              <wvf_param>{0}</wvf_param>\n'.format(smear_factor(y))+
                '              <wvfIntPar>{0}</wvfIntPar>\n'.format(smear_iter(y))+
                '            <no_smear_dir>3</no_smear_dir>\n'+
                '          </SmearingParam>\n'+
                '          <LinkSmearing>\n'+
                '            <LinkSmearingType>APE_SMEAR</LinkSmearingType>\n'+
                '            <link_smear_fact>{0}</link_smear_fact>\n'.format(smear_factor(y))+
                '            <link_smear_num>{0}</link_smear_num>\n'.format(smear_iter(y))+
                '            <no_smear_dir>3</no_smear_dir>\n'+
                '          </LinkSmearing>\n'+
                '        </Sink>\n'+
                '      </PropSink>\n'+
                '      <NamedObject>\n'+
                '        <gauge_id>default_gauge_field</gauge_id>\n'+
                '        <prop_ids>\n'+
                '          <elem>sh_{0}_prop</elem>\n'.format(x)+
                '          <elem>sh_{0}_prop</elem>\n'.format(x)+
                '        </prop_ids>\n'+
                '        <seqsource_id>sh_{0}_sh_{1}_seqsource</seqsource_id>\n'.format(x,y)+
                '      </NamedObject>\n'+
                '    </elem>\n\n')
        f.write('   <elem>\n'+
                '      <Name>PROPAGATOR</Name>\n'+
                '      <Frequency>1</Frequency>\n'+
                '      <Param>\n'+
                '        <version>10</version>\n'+
                '        <quarkSpinType>FULL</quarkSpinType>\n'+
                '        <obsvP>false</obsvP>\n'+
                '        <numRetries>1</numRetries>\n'+
                '        <FermionAction>\n'+
                '         <FermAct>WILSON</FermAct>\n'+
                '         <Mass>{0}</Mass>\n'.format(m)+
                '         <AnisoParam>\n'+
                '           <anisoP>true</anisoP>\n'+
                '           <t_dir>3</t_dir>\n'+
                '           <xi_0>2.0</xi_0>\n'+
                '           <nu>2.0</nu>\n'+
                '         </AnisoParam>\n'+
                '         <FermionBC>\n'+
                '           <FermBC>SIMPLE_FERMBC</FermBC>\n'+
                '           <boundary>1 1 1 -1</boundary>\n'+
                '         </FermionBC>\n'+
                '        </FermionAction>\n'+
                '    <InvertParam>\n'+
                '        <invType>QUDA_WILSON_INVERTER</invType>\n'+
                '        <WilsonParams>\n'+
                '            <Mass>{0}</Mass>\n'.format(m)+
                '            <AnisoParam>\n'+
                '              <anisoP>true</anisoP>\n'+
                '              <t_dir>3</t_dir>\n'+
                '              <xi_0>2.0</xi_0>\n'+
                '              <nu>2.0</nu>\n'+
                '            </AnisoParam>\n'+
                '         </WilsonParams>\n'+
                '          <RsdTarget>1e-8</RsdTarget>\n'+
                '          <RsdToleranceFactor>1000</RsdToleranceFactor>\n'+
                '          <Delta>0.5</Delta>\n'+
                '          <MaxIter>10000</MaxIter>\n'+
                '          <AntiPeriodicT>true</AntiPeriodicT>\n'+
                '          <SolverType>CG</SolverType>\n'+
                '          <Verbose>false</Verbose>\n'+
                '          <AsymmetricLinop>false</AsymmetricLinop>\n'+
                '          <CudaPrecision>SINGLE</CudaPrecision>\n'+
                '          <CudaReconstruct>RECONS_NONE</CudaReconstruct>\n'+
                '          <CudaSloppyPrecision>HALF</CudaSloppyPrecision>\n'+
                '          <CudaSloppyReconstruct>RECONS_NONE</CudaSloppyReconstruct>\n'+
                '          <AxialGaugeFix>false</AxialGaugeFix>\n'+
                '          <AutotuneDslash>true</AutotuneDslash>\n'+
                '        </InvertParam>\n'+
                '      </Param>\n'+
                '       <NamedObject>\n'+
                '        <gauge_id>default_gauge_field</gauge_id>\n'+
                '        <source_id>sh_{0}_sh_{1}_seqsource</source_id>\n'.format(x,y)+
                '        <prop_id>sh_{0}_sh_{1}_seqprop</prop_id>\n'.format(x,y)+
                '       </NamedObject>\n'+
                '    </elem>\n\n')
        f.write('   <elem>\n'+
                '      <Name>ERASE_NAMED_OBJECT</Name>\n'+
                '      <Frequency>1</Frequency>\n'+
                '      <NamedObject>\n'+
                '        <object_id>sh_{0}_sh_{1}_seqsource</object_id>\n'.format(x,y)+
                '      </NamedObject>\n'+
                '    </elem>\n\n')
        f.write(
                '    <elem>\n'+
                '      <Name>BAR3PTFN</Name>\n'+
                '      <Frequency>1</Frequency>\n'+
                '      <Param>\n'+
                '        <version>7</version>\n'+
                '        <j_decay>3</j_decay>\n'+
                '        <mom2_max>1</mom2_max>\n'+
                '      </Param>\n'+
                '      <NamedObject>\n'+
                '        <gauge_id>default_gauge_field</gauge_id>\n'+
                '        <prop_id>sh_{0}_prop</prop_id>\n'.format(x)+
                '        <bar3ptfn_file>bar3ptfn_sh_{0}_sh_{1}_CONFIGNAME.dat</bar3ptfn_file>\n'.format(x,y)+
                '        <seqprops>\n'+
                '          <elem>\n'+
                '            <seqprop_id>sh_{0}_sh_{1}_seqprop</seqprop_id>\n'.format(x,y)+
                '            <gamma_insertion>{0}</gamma_insertion>\n'.format(i)+
                '          </elem>\n'+
                '       </seqprops>\n'+
                '      </NamedObject>\n'+
                '    </elem>\n\n')
        f.write('   <elem>\n'+
                '      <Name>ERASE_NAMED_OBJECT</Name>\n'+
                '      <Frequency>1</Frequency>\n'+
                '      <NamedObject>\n'+
                '        <object_id>sh_{0}_sh_{1}_seqprop</object_id>\n'.format(x,y)+
                '      </NamedObject>\n'+
                '    </elem>\n\n')
    f.write('   <elem>\n'+
            '      <Name>ERASE_NAMED_OBJECT</Name>\n'+
            '      <Frequency>1</Frequency>\n'+
            '      <NamedObject>\n'+
            '        <object_id>sh_{0}_prop</object_id>\n'.format(x)+
            '      </NamedObject>\n'+
            '    </elem>\n\n')

f.write(
        '  </InlineMeasurements>\n'+
        '   <nrow>{0} {0} {0} {1}</nrow>\n'.format(s_size, t_size)+
        '</Param>\n'+
        '<RNG>\n'+
        '  <Seed>\n'+
        '    <elem>1</elem>\n'+
        '    <elem>2</elem>\n'+
        '    <elem>3</elem>\n'+
        '    <elem>4</elem>\n'+
        '  </Seed>\n'+
        '</RNG>\n'+
        '<Cfg>\n'+
        '  <cfg_type>SZINQIO</cfg_type>\n'+
        '  <cfg_file>CONFIGNAME.lime</cfg_file>\n'+
        '</Cfg>\n'+
        '</chroma>\n')